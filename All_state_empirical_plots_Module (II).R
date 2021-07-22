
#-----------Prortin DNA Map--------------------------------#

setwd("C:/protein/TDP43/main/")
set.seed(123)
library(tidyverse)
library(gplots)

system<-dir()[grep(".pdb", dir())]

filtered_inter_distance_list_DNA<-list()
filtered_inverse_inter_distance_list_DNA<-list()
distance_matrix<-list()
contact_matrix<-list()
protein2<-c("TSDLIVLGLPWKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRWCDCKLPNS")
protein_wt<-strsplit(protein2, "")[[1]]
DNA_seq<-c("G","T","T","G","A","G","C","G","T","T")

euclid<-function(list_item){
  
  A<-apply(list_item[,7:9], 2, as.numeric)
  
  #||A-B||^2 =  A^2 + B^2 - 2*A*B
  
  d <- sqrt(matrix(diag(A %*% t(A)), nrow=dim(A)[1], ncol=dim(A)[1])
            + t(matrix(diag(A %*% t(A)), nrow=dim(A)[1], ncol=dim(A)[1]))
            - 2*A %*% t(A) )
  return(d)
} #euclid distance function

contact<-function(list_item){
  
  A<-apply(list_item[,7:9], 2, as.numeric)
  
  #||A-B||^2 =  A^2 + B^2 - 2*A*B
  
  d <- sqrt(matrix(diag(A %*% t(A)), nrow=dim(A)[1], ncol=dim(A)[1])
            + t(matrix(diag(A %*% t(A)), nrow=dim(A)[1], ncol=dim(A)[1]))
            - 2*A %*% t(A) )
  
  d1<-ifelse(d<=10,1,0)
  return(d1)
}
  
  
  
for(kk in 1:4){
  
  print(paste("start", kk, "at", Sys.time(), sep=" "))
  
  
  data<-read.table(system[kk], sep="\t")
  
  start<- grep("MODEL", data[,1])
  end<-grep("ENDMDL", data[,1])
  
  protein_snap<-list()
  DNA_snap<-list()
  DNA_CM_snap<-list()
  snap<-list()
  
  for(i in seq_along(start)){
    
    current_snap <- str_squish(data[(start[i]+1):(end[i]-2),1])
    
    prot_snap <- current_snap[grep("CA", current_snap)]
    dna_snap <- current_snap[grep(" B ", current_snap)]
    
    protein_snap <- as.data.frame(do.call(rbind, strsplit(prot_snap, " ")))
    
    dna <- as.data.frame(do.call(rbind, strsplit(dna_snap, " ")))
    
    suppressWarnings(
      aa<-matrix(dna[,6], ncol=2,nrow=(nrow(dna)+1)))
    nucl_id<-dna[aa[,1] != aa[,2],4]
    
    snap_temp<-cbind(rep("ATOM", length(unique(dna$V6))),
                     rep("DNA", length(unique(dna$V6))),
                     rep("XX", length(unique(dna$V6))),
                     c(nucl_id),
                     rep("B", length(unique(dna$V6))),
                     (c(1:length(unique(dna$V6)))+78),
                     rep("x", length(unique(dna$V6))),
                     rep("y", length(unique(dna$V6))),
                     rep("z", length(unique(dna$V6))),
                     rep(1.00, length(unique(dna$V6))),
                     rep(0.00, length(unique(dna$V6))),
                     rep("H", length(unique(dna$V6))))
    
    
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
    snap_temp[,7:9]<-center_of_mass
    DNA_CM_snap <-as.data.frame(snap_temp)
    
    
    snap[[i]]<-rbind(protein_snap, DNA_CM_snap)
    
  }
  print("--------------------DNA and Protein Merged Successfully------------------------")
  
  #----------------------Distance--------------------#
  distance_matrix[[kk]] <-lapply(snap, euclid)
  contact_matrix[[kk]] <-lapply(snap, contact)
}


distance.matrix<-list()
contact.matrix<-list()
for(i in 1:4){ distance.matrix[[i]] <- Reduce("+", distance_matrix[[i]])/length(distance_matrix[[i]])
contact.matrix[[i]] <- Reduce("+", contact_matrix[[i]])/length(contact_matrix[[i]])}
#for(i in 1:4){ write.table(distance.matrix[[i]], paste0("distance_", system[i], ".txt", collapse = "")) }
names(distance.matrix)<- system
names(contact.matrix)<- system



#---------------------PLOTS---------------------------------#
#------------------------------------------------------------#
#------------------------------------------------------------#



#---------------Protein-Protein Average Inter-residual Distance-----------------------------------#


my_palette <-colorRampPalette(c("white","skyblue","blue"
                                ,"orange","red","black"))(n = 400)


for(i in 1:4){
  
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
  
  
}


#-----------------Protein-DNA Average Inter-residual Distance--------------------

my_palette <-colorRampPalette(c("white","skyblue","blue"
                                ,"orange","red","black"))(n = 400)


for(i in 1:4){
  
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







#------------------Contact Protein-Protein--------------#

my_palette <-colorRampPalette(c("pink",
                                "red", "black"))(n = 50)


for(k in 1:4){
  
  current_map<-contact.matrix[[k]][1:78,1:78]
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
            key.xlab = "Protein residues contact intensity",
            key.ylab = "Density",
            col= my_palette,
            font.lab=9,
            xlab = "Protein_Residue",
            ylab = "Protein_Residue",
            main = system[k],
            labCol = protein_wt,
            labRow = protein_wt)
  
  
}


#------------------Contact Protein-DNA--------------#

my_palette <-colorRampPalette(c("pink",
                                "red", "black"))(n = 50)


for(k in 1:4){
  
  current_map<-contact.matrix[[k]][77:88,1:78]
  
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
            main = system[k],
            labCol = protein_wt,
            labRow = DNA_seq)
  
  
}
















