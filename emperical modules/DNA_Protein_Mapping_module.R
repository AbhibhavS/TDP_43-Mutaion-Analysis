
#-----------Prortin DNA Map--------------------------------#

setwd("D:/Abhibhav/TDP43")
library(tidyverse)

files<-dir()[grep(".pdb", dir())]
my_map<-list()

for(file in seq_along(files)){
  
  data<-read.delim(files[file]  , sep="\t")
  
  start<- grep("MODEL", data[,1])
  end<-grep("ENDMDL", data[,1])
  
  protein_snap<-list()
  dna_snap<-list()

  for(i in seq_along(start)){
    
    current_snap <- str_squish(data[(start[i]+1):(end[i]-2),1])
    suppressWarnings(
      current_snap <- do.call(rbind, strsplit(current_snap, " ")))
    protein_snap[[i]] <- as.data.frame(current_snap[which(current_snap[,3]=="CA"),])
    dna_snap[[i]] <- as.data.frame(current_snap[grep("D", current_snap[,4]),])
    
    
  }
  
  
  center_of_mass_in_each_snap<-list()
  
  for(k in seq_along(dna_snap)){
    
    dna<-dna_snap[[k]]
    
    center_of_mass<-NULL
    for(j in seq_along(unique(dna$V6))){
      
      neucleotide <- dna[which(dna$V6 == unique(dna$V6)[j]), ]
      
      neucleotide_mass<-c(ifelse(c(1:nrow(neucleotide))%in%grep("H", neucleotide[,3]), 1, 0) +
                            ifelse(c(1:nrow(neucleotide))%in%grep("N", neucleotide[,3]), 14, 0) +
                            ifelse(c(1:nrow(neucleotide))%in%grep("P", neucleotide[,3]), 31, 0) +
                            ifelse(c(1:nrow(neucleotide))%in%grep("C", neucleotide[,3]), 12, 0) +
                            ifelse(c(1:nrow(neucleotide))%in%grep("O", neucleotide[,3]), 16, 0))
      
      cm_x<- sum(neucleotide_mass*as.numeric(neucleotide$V7))/sum(neucleotide_mass)
      cm_y<- sum(neucleotide_mass*as.numeric(neucleotide$V8))/sum(neucleotide_mass)
      cm_z<- sum(neucleotide_mass*as.numeric(neucleotide$V9))/sum(neucleotide_mass)
      
      center_of_mass<-rbind(center_of_mass, c(cm_x, cm_y, cm_z))
    }
    
    center_of_mass_in_each_snap[[k]]<-center_of_mass
  }
  
  
  map<-matrix(0, ncol=nrow(protein_snap[[1]]), 
              nrow=nrow(center_of_mass_in_each_snap[[1]]))
  
  
  for(k in seq_along(protein_snap)){
    
    for(p in 1:nrow(protein_snap[[k]])){
      for(d in 1:nrow(center_of_mass_in_each_snap[[k]])){
        
        #map[d,p]<- ifelse(sqrt(sum((as.numeric(protein_snap[[k]][p,7:9])-center_of_mass_in_each_snap[[k]][d,])^2)) <= 10, 
                          #map[d,p]+1, map[d,p]+0)
        map[d,p] <- map[d,p]+ sqrt(sum((as.numeric(protein_snap[[k]][p,7:9])-center_of_mass_in_each_snap[[k]][d,])^2))
        
      }
    }
  }
  
  my_map[[file]]<-map
}

#------------------------PLOT--------------------------------#

library(gplots)

for(k in seq_along(files)){
  
  map<- my_map[[k]]/(length(protein_snap))
  aa<-matrix(dna_snap[[1]][,6], ncol=2, 
             nrow=(length(dna_snap[[1]][,6])+1))
  aa<-dna_snap[[1]][aa[,1] != aa[,2],4]
  
  my_palette <-colorRampPalette(c("darkblue", "white"))(n = 50)
  heatmap.2(map, trace="none", dendrogram = "none",
            margins = c(5,9), 
            density.info=c("none"),
            key = T, Rowv = F, Colv = F,
            cexRow = 0.8, cexCol = 0.7, 
            key.xlab = "Avg. Distance (Angstrom)",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "DNA_Residue",
            main = files[k],
            labRow = paste(aa,1:length(aa)),
            labCol = paste(protein_snap[[1]][,4], 
                           1:nrow(protein_snap[[1]]))
  )

}
  

  
  
  
  






