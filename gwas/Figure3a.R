setwd("I:/Quinoa/Results/SNP_Calling_V2/Final_data_analysis/Rcodes_figures")
library(qqman)


DTF<-read.table("Cq2A_MAT_plots_emmax_P_K_DTF_all_blue_NO_dup.phe.ps.qqman",sep="\t",header=F)
DTM<-read.table("Cq2A_MAT_plots_emmax_P_K_DTM_all_blue_NO_dup.phe.ps.qqman",sep="\t",header=F)
PH<-read.table("Cq2A_MAT_plots_emmax_P_K_PH_all_blue_NO_dup.phe.ps.qqman",sep="\t",header=F)
PL<-read.table("Cq2A_MAT_plots_emmax_P_K_PL_all_blue_NO_dup.phe.ps.qqman",sep="\t",header=F)

colnames(DTF)<-c("Marker","Chr","Pos","p")
colnames(DTM)<-c("Marker","Chr","Pos","p")
colnames(PH)<-c("Marker","Chr","Pos","p")
colnames(PL)<-c("Marker","Chr","Pos","p")

png("Cq2A_DTF_all_blue_NO_dup.phe.ps.qqman.png",type="cairo",width=889,height=535,units="px")
manhattan(DTF,chr="Chr",bp="Pos", snp="Marker", p="p",ylim=c(0,10), main = "DTF (Cq2A)",
          suggestiveline = -log10(8.98E-7),genomewideline = -log10(1.66877489e-8),
          cex = 1.5,cex.axis =1.3,col = c("black", "#CC3333"))
dev.off()

#-------------------------------------------------------------------------------------------#

png("Cq2A_DTM_all_blue_NO_dup.phe.ps.qqman.png",type="cairo",width=889,height=535,units="px")
manhattan(DTM,chr="Chr",bp="Pos", snp="Marker", p="p",ylim=c(0,12), main = "DTM (Cq2A)",
          suggestiveline = -log10(8.98E-7),genomewideline = -log10(1.66877489e-8),
          cex = 1.5,cex.axis =1.3,col = c("black", "#CC3333"))
dev.off()

#-------------------------------------------------------------------------------------------#

png("Cq2A_PH_all_blue_NO_dup.phe.ps.qqman.png",type="cairo",width=889,height=535,units="px")
manhattan(PH,chr="Chr",bp="Pos", snp="Marker", p="p",ylim=c(0,10), main = "DTF (Cq2A)",
          suggestiveline = -log10(8.98E-7),genomewideline = -log10(1.66877489e-8),
          cex = 1.5,cex.axis =1.3,col = c("black", "#CC3333"))
dev.off()

#-------------------------------------------------------------------------------------------#

png("Cq2A_PL_all_blue_NO_dup.phe.ps.qqman.png",type="cairo",width=889,height=535,units="px")
manhattan(PL,chr="Chr",bp="Pos", snp="Marker", p="p",ylim=c(0,10), main = "DTF (Cq2A)",
          suggestiveline = -log10(8.98E-7),genomewideline = -log10(1.66877489e-8),
          cex = 1.5,cex.axis =1.3,col = c("black", "#CC3333"))
dev.off()