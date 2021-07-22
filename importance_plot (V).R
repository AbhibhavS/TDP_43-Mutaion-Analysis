
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

min_max<- function(x){
  (x-min(x))/(max(x)-min(x))
}

system_importance<-list()

for(kk in seq_along(systems)){
  
  Enet<-read.csv(system_mat[,kk][1], row.names = 1)
  RF<-read.csv(system_mat[,kk][4], row.names = 1)
  #pca<-read.csv(system_mat[,kk][3], row.names = 1)
  
  model<-list(Enet, RF)#, pca)
  important<-list()
  
  for(k in seq_along(model)){
    
    residue_number<-as.numeric(unlist(regmatches( rownames(model[[k]]), 
                                                  gregexpr("[[:digit:]]+",  
                                                           rownames(model[[k]])))))
    
    imp_inter<-cbind(model[[k]], matrix(residue_number, ncol=2, byrow = T))
    
    imp_residue<-matrix(1:length(sort(unique(residue_number))), 
                        nrow=length(unique(residue_number)),2)
    
    for(i in 1:nrow(imp_residue)){
      imp_residue[i,2]<-sum(imp_inter[which(imp_inter[,2:3]==i, arr.ind = T)[,1],1])
    }
    
    imp_residue[,2]<- min_max(imp_residue[,2])
    
    important[[k]]<-imp_residue[,2]
    
  }
  
  names(important)<-c("Enet", "RF")#, "PCA")
  
  system_importance[[kk]]<-important
  
  
  my.importance<- Reduce("+", important)/length(model)
  my.importance<-min_max(my.importance[1:78])
  maximums<-c(1:78)[ggpmisc:::find_peaks(my.importance[1:78])]
  top_10<-maximums[order(my.importance[maximums], decreasing = T)][1:15]
  
  
  
  
  plot(NA,pch=15, cex=2, ylim=c(0,1.1),  yaxt="n",
       xlim=c(1,78),xaxt="n",
       xlab="", ylab="")
  segments(top_10, rep(-1, length(top_10)), col="grey",
           top_10, rep(1, length(top_10)), lwd=2, lty=3)
  #abline(h=mean(my.importance[1:78]))
  par(new=T)
  
  
  plot(1:78, my.importance[1:78], ylim=c(0,1.1), xlim=c(1,78),
       type="l", col="Blue", lwd=5, xaxt="n", yaxt="n",
       ylab="Importance", xlab="residues", main=paste(systems[kk],
                                                      "Ave. Imp/res (Enet+RF)", sep=" "))
  axis(1, at=1:78, labels = 1:78, 
       font = 2, tick = T, cex.axis=0.7)
  axis(2, at=seq(0,1,0.1), labels = seq(0,1,0.1), font = 2, tick = T,
       las=1)
  box(lwd=2)
  par(new=T)
  
  color<-rep("black", length(protein2))
  color[mut[,kk]]<-"Red"
  text(x=1:78, y=rep(1.1,78), 
       labels = protein2, 
       cex=0.8, font=9, col=color)
  par(new=T)
  plot(sheet, y=rep(1.05, length(sheet)), 
       pch=15, cex=2, ylim=c(0,1.1), yaxt="n",
       xlim=c(1,78), col="yellow", xaxt="n",
       xlab="", ylab="")
  par(new=T)
  plot(helix, y=rep(1.05, length(helix)), 
       pch=15, cex=2, ylim=c(0,1.1),  yaxt="n",
       xlim=c(1,78), col="red", xaxt="n",
       xlab="", ylab="")
}







#--------- only use to plot for separate model ---------------#



for(i in 2:1){
  
  plot(1:78, system_importance[[1]][[i]][1:78], ylim=c(0,1.1), xlim=c(1,78),
       type="l", col=colo[i], lwd=5, xaxt="n", yaxt="n",
       ylab="Importance", xlab="residues", main="Residue_importance")
  axis(1, at=1:78, labels = 1:78, 
       font = 2, tick = T, cex=4)
  axis(2, at=seq(0,1,0.1), labels = seq(0,1,0.1), font = 2, tick = T,
       las=1)
  
  box(lwd=2)
  par(new=T)
  
}

sheet<-c(4:7,28:35, 42:50, 64:66, 69:74)
helix<-c(14:23, 52:61)

text(x=1:78, y=rep(1.1,78), 
     labels = protein2, 
     cex=0.8, font=9)
par(new=T)
plot(sheet, y=rep(1.05, length(sheet)), 
     pch=15, cex=2, ylim=c(0,1.1), yaxt="n",
     xlim=c(1,78), col="yellow", xaxt="n",
     xlab="", ylab="")
par(new=T)
plot(helix, y=rep(1.05, length(helix)), 
     pch=15, cex=2, ylim=c(0,1.1),  yaxt="n",
     xlim=c(1,78), col="red", xaxt="n",
     xlab="", ylab="")



names(system_importance)<- systems

global<-cbind(c(protein2,DNA),
             system_importance$Global$Enet, 
             system_importance$Global$RF,
             Reduce("+",system_importance$Global)/2)


D169G<-cbind(c(protein2,DNA),
             system_importance$D169G$Enet, 
             system_importance$D169G$RF,
             Reduce("+",system_importance$D169G)/2)

D169G_I168A<-cbind(c(protein2,DNA),
             system_importance$`D169G-I168A`$Enet, 
             system_importance$`D169G-I168A`$RF,
             Reduce("+",system_importance$`D169G-I168A`)/2)

I168A<-cbind(c(protein2,DNA),
             system_importance$I168A$Enet, 
             system_importance$I168A$RF,
             Reduce("+",system_importance$I168A)/2)

D169G<-as.data.frame(D169G)
D169G_I168A<-as.data.frame(D169G_I168A)
I168A<-as.data.frame(I168A)
global<-as.data.frame(global)

names(D169G)<-c("atom", "Enet", "RF", "Average")
names(D169G_I168A)<-c("atom", "Enet", "RF", "Average")
names(I168A)<-c("atom", "Enet", "RF", "Average")
names(global)<-c("atom", "Enet", "RF", "Average")

#write.table(D169G, "Per_Residue_IMP_D169G.txt", row.names = F)
#write.table(D169G_I168A, "Per_Residue_IMP_D169G-I168A.txt", row.names = F)
#write.table(I168A, "Per_Residue_IMP_I168A.txt", row.names = F)
#write.table(global, "Per_Residue_IMP_global.txt", row.names = F)








