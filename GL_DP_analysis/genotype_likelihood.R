setwd("C:/Users/dilan/OneDrive/eLife/new revisions 3.1.2022/new_analyisis/genotyope_likelyhood/")

library(tidyverse)
library(patchwork)
library(ggpubr)


S2<-read.table("S2.txt",header=T,sep="\t")

GQ<-read.table("Hetero_GQ_mean.txt",header=F,sep="\t")

head(S2)

colnames(GQ)<-"GQ_Het"

s2_final<- S2 %>% mutate(het_persentage=(HET_SNPs/(HET_SNPs+REF_SNPs+ALT_SNPs)*100),
              meanDP_filterd=(REF_DP+ALT_DP+HET_DP)/3)

dim(S2)
dim(GQ)

GQ_2<-cbind(s2_final,GQ)

#write.table(s2_final,"s2_final.txt",row.names = F,quote = F,sep = "\t")

head(GQ_2)

#-----------------------------------------------------------------------#


a<-ggscatter(GQ_2, x = "DP_unfiltered", y = "het_persentage",
             add = "reg.line", conf.int = TRUE,
             cor.coef = T, cor.method = "pearson",
             xlab = "mean DP unfiltered", ylab = "heterozygosity",
             title = "")+
  theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))
  

b<-ggscatter(GQ_2, x = "HET_DP", y = "het_persentage",
             add = "reg.line", conf.int = TRUE,
             cor.coef = T, cor.method = "pearson",
             xlab = "mean DP HET", ylab = "heterozygosity",
             title = "")+
  theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))


c<-ggscatter(s2_final, x = "ALT_DP", y = "het_persentage",
             add = "reg.line", conf.int = TRUE,
             cor.coef = T, cor.method = "pearson",
             xlab = "mean DP ALT_HOMO", ylab = "heterozygosity",
             title = "")+
  theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

d<-ggscatter(s2_final, x = "REF_DP", y = "het_persentage",
             add = "reg.line", conf.int = TRUE,
             cor.coef = T, cor.method = "pearson",
             xlab = "mean DP REF_HOMO", ylab = "heterozygosity",
             title = "")+
  theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))

e<-ggscatter(s2_final, x = "meanDP_filterd", y = "het_persentage",
            add = "reg.line", conf.int = TRUE,
            cor.coef = T, cor.method = "pearson",
            xlab = "mean DP filtered", ylab = "heterozygosity",
            title = "")+
  theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))



a+b+c+d+e

#-----------------------------------------------------------------------#

ggscatter(GQ_2, x = "GQ_Het", y = "het_persentage",
             add = "reg.line", conf.int = TRUE,
             cor.coef = T, cor.method = "pearson",
             xlab = "mean GQ Het", ylab = "heterozygosity",
             title = "")+
  theme(text = element_text(size=14),plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))




#-----------------------------------------------------------------------#


S2_less5<-filter(S2,DP_unfiltered<5)

S2_high5<-filter(S2,DP_unfiltered>5)

dim(S2_less5)
dim(S2_high5)

par(mfrow = c(3, 2))

hist(S2$DP_unfiltered,main = "Histogram of mean DP (unfiltered SNPs)",
          xlab = "mean DP \n (unfiltered SNPs)",
     ylab = "Frequency")

median(S2$DP_unfiltered)

boxplot(S2_high5$GQ_Homo,S2_less5$GQ_Homo, names = c("DP>5","DP<5"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples (unfiltered) \n homozygous SNPs",cex.main=1)


boxplot(S2_high5$Hetero_GQ_mean,S2_less5$Hetero_GQ_mean,
        names = c("DP>5","DP<5"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples (unfiltered) \n heterozygous SNPs",cex.main=1)

#-----------------------------------------------------------------------#
mean(S2$meanDP_filterd)

mean(S2$Hetero_GQ_mean)

mean(S2$GQ_Homo)

S2_less5_filtered<-filter(S2,meanDP_filterd<8.4)

S2_high5_filtered<-filter(S2,meanDP_filterd>8.4)

hist(S2$meanDP_filterd,main = "Histogram of mean DP (filtered SNPs)",
     xlab = "mean DP \n (filtered SNPs)",
     ylab = "Frequency")


boxplot(S2_high5_filtered$GQ_Homo,S2_less5_filtered$GQ_Homo, 
        names = c("DP>8.4","DP<8.4"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples (filtered) \n homozygous SNPs",cex.main=1)

boxplot(S2_high5_filtered$Hetero_GQ_mean,S2_less5_filtered$Hetero_GQ_mean,
        names = c("DP>8.6","DP<8.6"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples (filtered) \n heterozygous SNPs",cex.main=1)


#-----------------------------------------------------------------------#

par(mfrow = c(2, 2))

boxplot(S2_high5$ALT_Homozygous_GQ_mean,S2_less5$ALT_Homozygous_GQ_mean, 
        names = c("DP>5","DP<5"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples \n (unfiltered) ALT homozygous SNPs",cex.main=1)


boxplot(S2_high5$REF_Homo_GQ,S2_less5$REF_Homo_GQ,
        names = c("DP>5","DP<5"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples\n (unfiltered) REF homozygous SNPs",cex.main=1)

boxplot(S2_high5_filtered$ALT_Homozygous_GQ_mean ,S2_less5_filtered$ALT_Homozygous_GQ_mean , 
        names = c("DP>8.4","DP<8.4"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples \n (filtered) ALT homozygous SNPs",cex.main=1)


boxplot(S2_high5_filtered$REF_Homo_GQ,S2_less5_filtered$REF_Homo_GQ, 
        names = c("DP>8.4","DP<8.4"),
        ylab="GQ",
        main="Comparison of GQ between high and low DP samples \n (filtered) REF homozygous SNPs",cex.main=1)


