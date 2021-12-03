setwd("I:/Quinoa/Results/SNP_Calling_V2/Final_data_analysis/PCA/seperate_chromsomes")

PCA<-read.csv("combined_PC_sep_chromosomes.csv")

head(PCA)
dim(PCA)

id<-read.csv("../pca_all_no_dup_GWAS.csv")

id<-id[,1]

cbind(id,PCA)->PCA

pacman::p_load(ggplot2,tidyverse,patchwork) 

(a+b+c)/(d+e+f)/(g+h+i)/(j+k+l)/(m+n+o)/(p+q+r)

#------------------------------------------------------------------------------------------------#
a<-ggplot(PCA,aes(PC1_1,PC2_1))+geom_point()+
  scale_x_continuous(name="PC1 (25.33%)")+
  scale_y_continuous(name="PC2 (10.61%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq1A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
b<-ggplot(PCA,aes(PC1_2,PC2_2))+geom_point()+
  scale_x_continuous(name="PC1 (23.47%)")+
  scale_y_continuous(name="PC2 (9.29%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq2A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
c<-ggplot(PCA,aes(PC1_3,PC2_3))+geom_point()+
  scale_x_continuous(name="PC1 (20.87%)")+
  scale_y_continuous(name="PC2 (8.43%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq3A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
d<-ggplot(PCA,aes(PC1_4,PC2_4))+geom_point()+
  scale_x_continuous(name="PC1 (20.56%)")+
  scale_y_continuous(name="PC2 (9.77%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq4A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
e<-ggplot(PCA,aes(PC1_5,PC2_5))+geom_point()+
  scale_x_continuous(name="PC1 (25.85%)")+
  scale_y_continuous(name="PC2 (11.01%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq5A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
f<-ggplot(PCA,aes(PC1_6,PC2_6))+geom_point()+
  scale_x_continuous(name="PC1 (22.20%)")+
  scale_y_continuous(name="PC2 (13.60%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq6A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
g<-ggplot(PCA,aes(PC1_7,PC2_7))+geom_point()+
  scale_x_continuous(name="PC1 (22.54%)")+
  scale_y_continuous(name="PC2 (9.59%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq7A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
h<-ggplot(PCA,aes(PC1_8,PC2_8))+geom_point()+
  scale_x_continuous(name="PC1 (26.11%)")+
  scale_y_continuous(name="PC2 (11.30%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq8A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
i<-ggplot(PCA,aes(PC1_9,PC2_9))+geom_point()+
  scale_x_continuous(name="PC1 (25.15%)")+
  scale_y_continuous(name="PC2 (9.61%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq9A")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
j<-ggplot(PCA,aes(PC1_10,PC2_10))+geom_point()+
  scale_x_continuous(name="PC1 (22.37%)")+
  scale_y_continuous(name="PC2 (10.77%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq1B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
k<-ggplot(PCA,aes(PC1_11,PC2_11))+geom_point()+
  scale_x_continuous(name="PC1 (24.31%)")+
  scale_y_continuous(name="PC2 (9.08%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq2B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
l<-ggplot(PCA,aes(PC1_12,PC2_12))+geom_point()+
  scale_x_continuous(name="PC1 (26.94%)")+
  scale_y_continuous(name="PC2 (12.83%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq3B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
m<-ggplot(PCA,aes(PC1_13,PC2_13))+geom_point()+
  scale_x_continuous(name="PC1 (23.83%)")+
  scale_y_continuous(name="PC2 (12.30%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq4B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
n<-ggplot(PCA,aes(PC1_14,PC2_14))+geom_point()+
  scale_x_continuous(name="PC1 (22.73%)")+
  scale_y_continuous(name="PC2 (12.78%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq5B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
o<-ggplot(PCA,aes(PC1_15,PC2_15))+geom_point()+
  scale_x_continuous(name="PC1 (28.55%)")+
  scale_y_continuous(name="PC2 (9.15%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq6B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
p<-ggplot(PCA,aes(PC1_16,PC2_16))+geom_point()+
  scale_x_continuous(name="PC1 (24.75%)")+
  scale_y_continuous(name="PC2 (8.21%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq7B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
q<-ggplot(PCA,aes(PC1_17,PC2_17))+geom_point()+
  scale_x_continuous(name="PC1 (26.35%)")+
  scale_y_continuous(name="PC2 (11.28%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq8B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#
r<-ggplot(PCA,aes(PC1_18,PC2_18))+geom_point()+
  scale_x_continuous(name="PC1 (27.19%)")+
  scale_y_continuous(name="PC2 (8.29%)")+
  theme_classic()+
  theme(text=element_text(family="Arial", size=12, color = "black"),legend.position = "none")+
  labs(title="Cq9B")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(color="Black", size=14, face="bold.italic"))+
  theme(axis.text.x = element_text(colour = "black", size = 11))+
  theme(axis.text.y = element_text(colour = "black", size = 11))+
  theme(axis.line = element_line(colour = "black", size = 0.60, linetype = "solid"))+
  theme(panel.background = element_rect(colour = "black", size = 0.60))

#------------------------------------------------------------------------------------------------#


