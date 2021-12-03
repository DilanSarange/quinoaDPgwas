
setwd("F:/data_analysis_from_the_server/summary_stat")

library(tidyverse)
library(CMplot)
library(ggplot2)
theme_set(theme_light()+ theme(legend.position = "right"))

#-------------------------------------------------------------------------------------------------------------------#

FST_01<-read.table("input/quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.maf_0.01.recode.vcf.gz_FST_higlandvslowland_10KB.windowed.weir.fst",sep = "\t",header = T,na.strings = c("na","nan","NA"))

head(FST_01)

FST_01 %>% mutate(
  SNP = paste("S",CHROM, BIN_END, sep = "_"))->FST_01


FST_01<-FST_01[,c(7,1,3,5,6)]

colnames(FST_01)<-c("SNP","chr","pos","Weighted_FST","Mean_FST")

quantile(FST_01$Weighted_FST,0.99,na.rm=T)

summary(FST_01)

#----------------------------------------------------------------------------------------------------------------#

Tajima_01_high<-read.table("input/TajimaD_window_10KBhighland.pop_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.maf_0.01.recode.vcf.Tajima.D",sep = "\t",header = T,na.strings = c("na","nan","NA"))

head(Tajima_01_high)
summary(Tajima_01_high)

mean(Tajima_01_high$TajimaD,na.rm = T)

Tajima_01_high %>% mutate(
  SNP = paste("S",CHROM, BIN_START, sep = "_"))-> Tajima_01_high

Tajima_01_high<-Tajima_01_high[,c(5,1,2,4)]

colnames(Tajima_01_high)<-c("SNP","chr","pos","TajimaD_high")

summary(Tajima_01_high)

#----------------------------------------------------------------------------------------------------------------#

Tajima_01_low<-read.table("input/TajimaD_window_10KBlowland.pop_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.maf_0.01.recode.vcf.Tajima.D",sep = "\t",header = T,na.strings = c("na","nan","NA"))

head(Tajima_01_low)

summary(Tajima_01_low)

mean(Tajima_01_low$TajimaD,na.rm = T)

Tajima_01_low %>% mutate(
  SNP = paste("S",CHROM, BIN_START, sep = "_"))-> Tajima_01_low

Tajima_01_low<-Tajima_01_low[,c(5,1,2,4)]

colnames(Tajima_01_low)<-c("SNP","chr","pos","TajimaD_low")

summary(Tajima_01_low)

#----------------------------------------------------------------------------------------------------------------#

Pi_01_high<-read.table("input/Pi_window_10KBhighland.pop_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.maf_0.01.recode.vcf.windowed.pi",sep = "\t",header = T,na.strings = c("na","nan","NA"))

head(Pi_01_high)

summary(Pi_01_high)

Pi_01_high %>% mutate(
  SNP = paste("S",CHROM, BIN_END, sep = "_"))-> Pi_01_high

Pi_01_high<-Pi_01_high[,c(6,1,3,5)]

colnames(Pi_01_high)<-c("SNP","chr","pos","Pi_high")

#----------------------------------------------------------------------------------------------------------------#

Pi_01_low<-read.table("input/Pi_window_10KBlowland.pop_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.maf_0.01.recode.vcf.windowed.pi",sep = "\t",header = T,na.strings = c("na","nan","NA"))

head(Pi_01_low)

summary(Pi_01_low)

Pi_01_low %>% mutate(
  SNP = paste("S",CHROM, BIN_END, sep = "_"))-> Pi_01_low

Pi_01_low<-Pi_01_low[,c(6,1,3,5)]

colnames(Pi_01_low)<-c("SNP","chr","pos","Pi_low")

#----------------------------------------------------------------------------------------------------------------#

left_join(Pi_01_high,Pi_01_low,by=c("SNP","chr","pos")) %>% 
  left_join(.,Tajima_01_high,by=c("SNP","chr","pos")) %>% 
    left_join(.,Tajima_01_low,by=c("SNP","chr","pos")) %>% 
      left_join(.,FST_01,by=c("SNP","chr","pos")) -> SNP_summary_stat

SNP_summary_stat %>% mutate(Pi_high_ratio_low=Pi_high/Pi_low) %>% mutate(Pi_low_ratio_high=Pi_low/Pi_high)->SNP_summary_stat

#write.csv(SNP_summary_stat,"SNP_summary_stat.csv",row.names = F,quote = F)

SNP_summary_stat<-read.csv("SNP_summary_stat.csv",na.strings = "NA")
#----------------------------------------------------------------------------------------------------------------#

head(SNP_summary_stat)

quantile(SNP_summary_stat$Pi_high,0.99,na.rm=T)

CMplot(SNP_summary_stat[,c(1,2,3,4)], type="h",plot.type="m", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylab=expression(italic(pi)[" (Type I)"]),ylim=NULL,outward=TRUE,threshold=0.002962756 ,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

#----------------------------------------------------------------------------------------------------------------#

head(SNP_summary_stat)

quantile(SNP_summary_stat$Pi_low,0.99,na.rm=T)

CMplot(SNP_summary_stat[,c(1,2,3,5)], type="h",plot.type="m", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylab=expression(italic(pi)[" (Type II)"]),ylim=NULL,outward=TRUE,threshold=0.003050315 ,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

#----------------------------------------------------------------------------------------------------------------#

head(SNP_summary_stat)

quantile(SNP_summary_stat$Pi_high_ratio_low,0.99,na.rm=T)

CMplot(SNP_summary_stat[,c(1,2,3,11)], type="h",plot.type="m", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylab=expression(italic(pi)[" (TypeI / Type II)"]),ylim=NULL,outward=TRUE,threshold=5.257905 ,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

#----------------------------------------------------------------------------------------------------------------#

head(SNP_summary_stat)

quantile(SNP_summary_stat$Pi_low_ratio_high,0.99,na.rm=T)

CMplot(SNP_summary_stat[,c(1,2,3,12)], type="h",plot.type="m", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylab=expression(italic(pi)[" (TypeII / Type I)"]),ylim=NULL,outward=TRUE,threshold=4.35783 ,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


#----------------------------------------------------------------------------------------------------------------#

head(SNP_summary_stat)

quantile(SNP_summary_stat$Weighted_FST,0.99,na.rm=T)

CMplot(SNP_summary_stat[,c(1,2,3,8)], type="h",plot.type="m", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylab=expression(italic("Fst")["(Type I vs Type II)"]),ylim=NULL,outward=TRUE,threshold=0.8085537 ,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

#----------------------------------------------------------------------------------------------------------------#

head(SNP_summary_stat)

summary(SNP_summary_stat)

quantile(SNP_summary_stat$Mean_FST,0.99,na.rm=T)

CMplot(SNP_summary_stat[,c(1,2,3,9)], type="h",plot.type="m", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylab=expression(italic("Fst")["(Type I vs Type II)"]),ylim=NULL,outward=TRUE,threshold=0.6239069 
 ,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


#----------------------------------------------------------------------------------------------------------------#

# plottting data 

CMplot(SNP_summary_stat, type="h",plot.type="c", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylim=NULL,outward=TRUE,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

CMplot(SNP_summary_stat[c(1,2,3,6,7,4,5,8,11)], type="h",plot.type="c", 
       chr.labels=paste("Cq",c("1A","2A","3A","4A","5A","6A","7A","8A","9A","1B","2B","3B","4B","5B","6B","7B","8B","9B"), sep=""),band=0.5, LOG10=F, 
       ylim=NULL,outward=TRUE,threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.4,chr.den.col=NULL,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


summary(SNP_summary_stat)


t.test(SNP_summary_stat$Pi_high,SNP_summary_stat$Pi_low)
summary(SNP_summary_stat)

mean(SNP_summary_stat$Pi_high,na.rm = T)/mean(SNP_summary_stat$Pi_low,na.rm = T)

library(ggridges)
library(patchwork)

gather(SNP_summary_stat,"Pi_high","Pi_low",key="SNP_summary",value = "value")->X

head(X)

X %>% 
ggplot(aes(y=SNP_summary,x=value,fill=SNP_summary))+
  geom_density_ridges(scale=0.8)


X %>% 
  ggplot(aes(y=SNP_summary,x=value,fill=factor(stat(quantile))))+
  stat_density_ridges(scale=2,na.rm = T,geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      quantiles = 4, quantile_lines = TRUE
  ) +
  scale_fill_viridis_d(name = "Quantiles")  



gather(SNP_summary_stat,"TajimaD_high","TajimaD_low",key="SNP_summary2",value = "value2")->Y


Y %>% 
  ggplot(aes(y=SNP_summary2,x=value2,fill=factor(stat(quantile))))+
  stat_density_ridges(scale=2,na.rm = T,geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      quantiles = 4, quantile_lines = TRUE
  ) +
  scale_fill_viridis_d(name = "Quantiles")  
  


quantile(SNP_summary_stat$Pi_high,0.99,na.rm=T)

mean(SNP_summary_stat$Pi_high,na.rm=T)/ mean(SNP_summary_stat$Pi_low,na.rm=T)

a<-  ggplot(SNP_summary_stat,aes(Pi_high))+
  geom_density()+
  geom_vline(xintercept = 0.0006, linetype="dashed", 
             color = "red", size=1)+
  annotate(geom="text", x=.003, y=920, label="mean=0.00061",
           color="red")

mean(SNP_summary_stat$Pi_low,na.rm=T)

b<-  ggplot(SNP_summary_stat,aes(Pi_low))+
  geom_density()+
  geom_vline(xintercept = 0.00057, linetype="dashed", 
             color = "red", size=1)+
  annotate(geom="text", x=.003, y=920, label="mean=0.00057",
           color="red")

mean(SNP_summary_stat$TajimaD_high,na.rm=T)

c<-ggplot(SNP_summary_stat,aes(TajimaD_high))+
  geom_density()+
  geom_vline(xintercept = 1.02127, linetype="dashed", 
             color = "red", size=1)+
  annotate(geom="text", x=4, y=.20, label="mean=1.02127",
           color="red")

mean(SNP_summary_stat$TajimaD_low,na.rm=T)

d<-ggplot(SNP_summary_stat,aes(TajimaD_low))+
  geom_density()+
  geom_vline(xintercept = 0.7536, linetype="dashed", 
             color = "red", size=1)+
  annotate(geom="text", x=4, y=.2, label="mean=0.7536",
           color="red")


(a|c)/(b|d)


ls()
#-------------------------------------------------------------------------------------------------------------------------------#

#SNP_summary_stat<-read.csv("SNP_summary_stat.csv",na.strings = "NA")

str(summary_stat)

quantile(SNP_summary_stat$Weighted_FST,0.99,na.rm=T)
quantile(SNP_summary_stat$Mean_FST,0.99,na.rm=T)







