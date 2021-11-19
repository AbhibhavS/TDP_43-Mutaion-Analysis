
#-------------Functions-----------

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


mdml.matrix<-function(md.data, type){
  
  if(!type%in%c("MD","ML")){stop("Invalid Type: MD or ML")}
  
  start<- grep("MODEL", md.data[,1])
  end<-grep("ENDMDL", md.data[,1])
  
  protein_snap<-list()
  DNA_snap<-list()
  DNA_CM_snap<-list()
  snap<-list()
  
  print(paste("Total Snaps: ", as.character(length(start))))
  
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
  distance_matrix <-lapply(snap, euclid)
  
  if(type=="MD"){
    contact_matrix <-lapply(snap, contact)
    aa<-list(distance_matrix, contact_matrix)
    names(aa)<-c("distance_matrix", "contact_matrix")
    print("------list of distance and contact matrix Successfull--------")
    return(aa)
  }
  else if(type=="ML"){
    
    coord <- snap[[1]]
    interaction<-matrix(NA, nrow(coord), nrow(coord))
    
    for(i in 1:nrow(coord)){
      for(j in 1:nrow(coord)){
        
        interaction[i,j]<- paste0(coord$V4[i],"( atom-",as.character(i),"-)",
                                  "-", coord$V4[j],"(atom-",as.character(j),"-)", 
                                  collapse = "")
        
      }
    }
    
    interact <- interaction[upper.tri(interaction)]
    print("------interaction matrix Successfull--------")
    dist.mat<-lapply(distance_matrix, function(d){d[upper.tri(d)]})
    inter_distance<-do.call(rbind, dist.mat)
    colnames(inter_distance)<-interact
    #inverse_inter_distance<-1/inter_distance
    return(inter_distance)}
  
}



