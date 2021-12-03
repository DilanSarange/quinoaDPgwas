setwd("/Quinoa/Results/SNP_Calling_V2/Final_data_analysis/LD_decay")

library(ggplot2)

#-----------------------------------------------------------------------------------------------------#

A_genome<-read.table("LD_DECAY_plots_A_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.stat.gz.bin.gz",header = F,sep = "\t")

head(A_genome)

A_genome<-replace(A_genome,"Genome","A")

A_genome<-A_genome[,-3:-5]

colnames(A_genome)<-c("Bin","LD","No of pairs","Genome")
head(A_genome)

tail(A_genome)


#-----------------------------------------------------------------------------------------------------#

B_genome<-read.table("LD_DECAY_plots_B_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf0.05.stat.gz.bin.gz",header = F,sep = "\t")

B_genome<-replace(B_genome,"Genome","B")
head(B_genome)

B_genome<-B_genome[,-3:-5]

colnames(B_genome)<-c("Bin","LD","No of pairs","Genome")

both_genomes<-rbind(A_genome,B_genome)

tail(both_genomes)

View(B_genome)

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# plotting 

c<-ggplot() +geom_line(data = both_genomes,aes(x=Bin,y=LD,colour=Genome),size=1.5,alpha=1)+
  labs(x="Pairwise distance (kb)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,50*10^3,100*10^3,150*10^3,200*10^3,250*10^3,300*10^3),labels=c("0","50","100","150","200","250","300"))+
  theme_classic()+
  theme(legend.position=c(0.8,0.6),legend.text = element_text(colour = "black", size = 12),
        axis.text=element_text(size=12),axis.title=element_text(size=17))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))+
  scale_color_manual(name = "Genome", labels = c("A ", "B"),values=c('#003366','#CC0033'))+
  ylim(0,0.8)+
  geom_hline(yintercept=c(0.2,0.1), linetype="dashed", 
             color = "black", size=0.5)

#-----------------------------------------------------------------------------------------------------#

Cq1A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_01_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq1A<-replace(Cq1A,"Chromosome","Cq1A")

Cq1A<-Cq1A[,-3:-6]

colnames(Cq1A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq2A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_02_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq2A<-replace(Cq2A,"Chromosome","Cq2A")

Cq2A<-Cq2A[,-3:-6]

colnames(Cq2A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq3A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_03_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq3A<-replace(Cq3A,"Chromosome","Cq3A")

Cq3A<-Cq3A[,-3:-6]

colnames(Cq3A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq4A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_04_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq4A<-replace(Cq4A,"Chromosome","Cq4A")

Cq4A<-Cq4A[,-3:-6]

colnames(Cq4A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq5A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_05_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq5A<-replace(Cq5A,"Chromosome","Cq5A")

Cq5A<-Cq5A[,-3:-6]

colnames(Cq5A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq6A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_06_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq6A<-replace(Cq6A,"Chromosome","Cq6A")

Cq6A<-Cq6A[,-3:-6]

colnames(Cq6A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq7A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_07_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq7A<-replace(Cq7A,"Chromosome","Cq7A")

Cq7A<-Cq7A[,-3:-6]

colnames(Cq7A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq8A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_08_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq8A<-replace(Cq8A,"Chromosome","Cq8A")

Cq8A<-Cq8A[,-3:-6]

colnames(Cq8A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq9A<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_09_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq9A<-replace(Cq9A,"Chromosome","Cq9A")

Cq9A<-Cq9A[,-3:-6]

colnames(Cq9A)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq1B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_10_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq1B<-replace(Cq1B,"Chromosome","Cq1B")

Cq1B<-Cq1B[,-3:-6]

colnames(Cq1B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq2B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_11_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq2B<-replace(Cq2B,"Chromosome","Cq2B")

Cq2B<-Cq2B[,-3:-6]

colnames(Cq2B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq3B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_12_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq3B<-replace(Cq3B,"Chromosome","Cq3B")

Cq3B<-Cq3B[,-3:-6]

colnames(Cq3B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq4B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_13_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq4B<-replace(Cq4B,"Chromosome","Cq4B")

Cq4B<-Cq4B[,-3:-6]

colnames(Cq4B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq5B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_14_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq5B<-replace(Cq5B,"Chromosome","Cq5B")

Cq5B<-Cq5B[,-3:-6]

colnames(Cq5B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq6B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_15_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq6B<-replace(Cq6B,"Chromosome","Cq6B")

Cq6B<-Cq6B[,-3:-6]

colnames(Cq6B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq7B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_16_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq7B<-replace(Cq7B,"Chromosome","Cq7B")

Cq7B<-Cq7B[,-3:-6]

colnames(Cq7B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq8B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_17_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq8B<-replace(Cq8B,"Chromosome","Cq8B")

Cq8B<-Cq8B[,-3:-6]

colnames(Cq8B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#

Cq9B<-read.table("./Sep_chrom/LD_DECAY_plots_Cq_18_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

Cq9B<-replace(Cq9B,"Chromosome","Cq9B")

Cq9B<-Cq9B[,-3:-6]

colnames(Cq9B)<-c("Bin","LD","Chromosome")

#-----------------------------------------------------------------------------------------------------#
#binding rows

genome_A_Chr<-rbind(Cq1A,Cq2A,Cq3A,Cq4A,Cq5A,Cq6A,Cq7A,Cq8A,Cq9A)

a<-ggplot()+ geom_line(data = genome_A_Chr,aes(x=Bin,y=LD,colour=Chromosome),size=1,alpha=1)+
  labs(x="Pairwise distance (kb)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,50*10^3,100*10^3,150*10^3,200*10^3,250*10^3,300*10^3),labels=c("0","50","100","150","200","250","300"))+
  theme_classic()+
  theme(legend.position=c(0.8,0.63),legend.text = element_text(colour = "black", size = 10),
        axis.text=element_text(size=12),axis.title=element_text(size=14))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))+
   scale_color_manual(values=c('#003366','#CC0033',"#66CC00","#66CCFF","#FF9933","#FF9900","#330033","#CC33CC","#FF66CC"))+
  geom_hline(yintercept=0.2, linetype="dashed", 
             color = "black", size=0.5)+
  ylim(0,0.8)

genome_B_Chr<-rbind(Cq1B,Cq2B,Cq3B,Cq4B,Cq5B,Cq6B,Cq7B,Cq8B,Cq9B)

b<-ggplot()+ geom_line(data = genome_B_Chr,aes(x=Bin,y=LD,colour=Chromosome),size=1,alpha=1)+
  labs(x="Pairwise distance (kb)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,50*10^3,100*10^3,150*10^3,200*10^3,250*10^3,300*10^3),labels=c("0","50","100","150","200","250","300"))+
  theme_classic()+
  theme(legend.position=c(0.8,0.63),legend.text = element_text(colour = "black", size = 10),
        axis.text=element_text(size=12),axis.title=element_text(size=14))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))+
  scale_color_manual(values=c('#003366','#CC0033',"#66CC00","#66CCFF","#FF9933","#FF9900","#330033","#CC33CC","#FF66CC"))+
  geom_hline(yintercept=0.2, linetype="dashed", 
             color = "black", size=0.5)+
  ylim(0,0.8)


library(patchwork)

(a|b)/(c+plot_spacer())


#-----------------------------------------------------------------------------------------------------#
# lowland and highland clusters 


low<-read.table("LD_DECAY_plots_lowland_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

head(low)

low<-replace(low,"Population","Lowland")

low<-low[,-3:-6]

colnames(low)<-c("Bin","LD","Population")

#-----------------------------------------------------------------------------------------------------#

high<-read.table("LD_DECAY_plots_highland_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf_LD_decay.stat.gz.bin.gz",header = F,sep = "\t")

head(high)

high<-replace(high ,"Population","Highland")

high<-high[,-3:-6]

colnames(high)<-c("Bin","LD","Population")

#-----------------------------------------------------------------------------------------------------#

#wild<-read.table("LD_DECAY_plots_wild_chenopodium_snps_maxmiss0.5_minDP5_maf0.05.ANNOTATED.recode.vcf.gz_LD_decay.stat.gz.bin",header = F,sep = "\t")

#head()

#wild<-replace(wild ,"Population","Wild")

#wild<-wild[,-3:-6]

#colnames(wild)<-c("Bin","LD","Population")

#-----------------------------------------------------------------------------------------------------#
#binding rows

#all_populations<-rbind(low,high,wild)

High_low_populations<-rbind(low,high)

write.csv(High_low_populations,"LD_decay_seperate_pop.csv")




ggplot()+ geom_line(data = High_low_populations,aes(x=Bin,y=LD,colour=Population),size=2,alpha=1)+
  labs(x="Pairwise distance (kb)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,50*10^3,100*10^3,150*10^3,200*10^3,250*10^3,300*10^3),labels=c("0","50","100","150","200","250","300"))+theme_classic()+
  theme(legend.position=c(0.8,0.6),axis.text=element_text(size=12),axis.title=element_text(size=14))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))+
  scale_color_manual(values=c('#003366','#CC0033'))+
  geom_hline(yintercept=c(0.2,0.1), linetype="dashed", 
             color = "black", size=0.5)+
  ylim(0,0.8)
    
  



ggplot()+ geom_line(data = High_low_populations,aes(x=Bin,y=LD,colour=Population),size=1,alpha=1)+labs(x="Pairwise distance (kb)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,50*10^3,100*10^3,150*10^3,200*10^3,250*10^3),labels=c("0","50","100","150","200","250"))+theme_classic()+
  theme(legend.position=c(0.8,0.6),axis.text=element_text(size=12),axis.title=element_text(size=14))+
  scale_color_manual(values=c('#003366','#CC0033'))+
  geom_hline(yintercept=c(0.2,0.1), linetype="dashed", 
             color = "black", size=0.5)+
  ylim(0,0.8)




