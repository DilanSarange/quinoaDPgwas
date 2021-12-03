setwd("C:/Users/Dilan/Desktop/data_analysis_from_the_server/ADMIXTURE")
setwd("I:/Quinoa/Results/SNP_Calling_V2/Final_data_analysis/ADMIXTURE/greedy_CLUMPP")
# examples 
#https://luisdva.github.io/rstats/model-cluster-plots/ 
#https://stackoverflow.com/questions/31570248/how-do-you-create-a-stacked-barplot-with-x-labels-and-borders-grouped-by-a-facto 
#https://stackoverflow.com/questions/41679888/how-to-group-order-data-in-r-for-a-barplot 


CV_error<-read.csv("/CV_error_final.csv")
head (CV_error)

library(ggplot2)
library(tidyverse)



ggplot(data=CV_error,aes(K,Error))+geom_line()+geom_point()+labs(title="Cross Validation Error",x="K",y="Cross-validation error")+
  scale_x_continuous(breaks=c(2:10),limits=c(2,10))+
  theme(plot.title = element_text(hjust = 0.5))

#------------------------------------------------------------------------------------------------------------------------------#

fam<-read.delim("../hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minMEANDP5_0.01.fam",sep = " ",header = F)

sample_ID<-data.frame(fam[,1])
dim(sample_ID)

colnames(sample_ID)<-"id"

id<-read.csv("../id_and_accessions.csv")

# conutry3.csv is used to remove countries with single entries

country<-read.csv("I:/Quinoa/Results/SNP_Calling_V2/Final_data_analysis/PCA/country3.csv")
type<-read.csv("I:/Quinoa/Results/SNP_Calling_V2/Final_data_analysis/Gene_identification/PCA/type_of_pop.csv")
#rm(sample_ID)

left_join(sample_ID,country,by="id") %>% left_join(.,type,by="id") %>% left_join(.,id,by="id") -> sample_ID

head(sample_ID)

#write.csv(sample_ID,"sample_ID.csv")

sample_ID<-read.csv("sample_ID.csv")

#------------------------------------------------------------------------------------------------------------------------------#

K_2=read.table("Q02/Q02.outfile")

K_2<-K_2[,c(2,3)]

dim(K_2)

K_2<-cbind(K_2,sample_ID)

write.csv(plot_data_K_2,"plot_data_K_2.csv")


plot_data_K_2 <- K_2 %>% 
  gather('Population', 'Ancestry', V2:V3) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_2, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))
  
  
  
  facet_grid(. ~ country, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  
  
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
    theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))

#------------------------------------------------------------------------------------------------------------------------------#

K_3=read.table("Q03/Q03.outfile")

K_3<-K_3[,c(2,3,4)]

K_3<-cbind(K_3,sample_ID)


plot_data_K_3 <- K_3 %>% 
  gather('Population', 'Ancestry', V2:V4) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_3, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ country, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  
  
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))


#------------------------------------------------------------------------------------------------------------------------------#

K_4=read.table("Q04/Q04.outfile")

K_4<-K_4[,c(2,3,4,5)]

K_4<-cbind(K_4,sample_ID)



plot_data_K_4 <- K_4 %>% 
  gather('Population', 'Ancestry', V2:V5) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_4, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ country, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))


#------------------------------------------------------------------------------------------------------------------------------#

K_5=read.table("Q05/Q5.outfile")

K_5<-K_5[,c(2,3,4,5,6)]

K_5<-cbind(K_5,sample_ID)


plot_data_K_5 <- K_5 %>% 
  gather('Population', 'Ancestry', V2:V6) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_5, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ country, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))


#------------------------------------------------------------------------------------------------------------------------------#

K_6=read.table("Q06/Q06.outfile")

K_6<-K_6[,c(2,3,4,5,6,7)]

K_6<-cbind(K_6,sample_ID)

plot_data_K_6 <- K_6 %>% 
  gather('Population', 'Ancestry', V2:V7) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_6, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ country, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))



#------------------------------------------------------------------------------------------------------------------------------#

K_7=read.table("Q07/Q07.outfile")

K_7<-K_7[,c(2,3,4,5,6,7,8)]

K_7<-cbind(K_7,sample_ID)


plot_data_K_7 <- K_7%>% 
  gather('Population', 'Ancestry', V2:V8) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

  
ggplot(plot_data_K_7, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ country, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  theme_classic()+facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))



tbl_df(plot_data_K_7)

x_K7<-data.frame(plot_data_K_7)

pop_V7<-x_K7[,c(2,5)]

tbl_df(pop_V7)

pop_V7 %>% distinct(id, .keep_all = TRUE)-> pop_v7_likely_assignment

write.csv(pop_V7_genotypes,"pop_genotypes_amixture_Q7.csv",row.names = F)

# id merging 

id<-read.csv("C:/Dilan/SNPcallingV2/Sample info/id_comparison.csv")

head(id)

left_join(id,pop_v7_likely_assignment,by="id")-> pop_V7_genotypes

ggplot_build(P)$data

#------------------------------------------------------------------------------------------------------------------------------#

K_8=read.table("Q08/Q08.outfile")

K_8<-K_8[,c(2,3,4,5,6,7,8,9)]

K_8<-cbind(K_8,sample_ID)

write.csv(plot_data_K_8,"plot_data_K_8.csv")


plot_data_K_8 <- K_8 %>% 
  gather('Population', 'Ancestry', V2:V9) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_8, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ Type, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))


# filtering K8

wild_k8<-filter(K_8,Type=="wild")


plot_data_wild_K_8 <- wild_k8 %>% 
  gather('Population', 'Ancestry', V2:V9) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))



ggplot(plot_data_wild_K_8, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ Type, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))

# filtering K8

highland_k8<-filter(K_8,Type=="highland")


plot_data_highland_K_8 <- highland_k8 %>% 
  gather('Population', 'Ancestry', V2:V9) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))



ggplot(plot_data_highland_K_8, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  theme(axis.text.x=element_text(angle=90, hjust=1,size=5))
  
  facet_grid(. ~ Type, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  # filtering K8
  
  lowland_k8<-filter(K_8,Type=="lowland")
  
  
  plot_data_lowland_K_8 <- lowland_k8 %>% 
    gather('Population', 'Ancestry', V2:V9) %>% 
    group_by(id) %>% 
    mutate(likely_assignment = Population[which.max(Ancestry)],
           assingment_prob = max(Ancestry)) %>% 
    arrange(likely_assignment, desc(assingment_prob)) %>% 
    ungroup() %>% 
    mutate(id = forcats::fct_inorder(factor(id)))
  
  
  
  ggplot(plot_data_lowland_K_8, aes(id, Ancestry, fill = Population)) +
    geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
    theme(axis.text.x=element_text(angle=90, hjust=1,size=5))
  
  facet_grid(. ~ Type, drop=TRUE, space="free", scales="free") +
    theme(panel.grid=element_blank()) +
    theme(panel.background=element_rect(fill=NA, colour="grey25")) +
    theme(panel.margin.x=grid:::unit(0.5, "lines")) +
    theme(axis.title.x=element_blank()) +
    theme(axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    theme(strip.background=element_blank()) +
    theme(strip.text=element_text(size=12)) +
    theme(legend.key=element_rect(colour="grey25")) +
    scale_x_discrete(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  
  
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))


#------------------------------------------------------------------------------------------------------------------------------#

K_9=read.table("Q09/Q09.outfile")


K_9<-K_9[,c(2,3,4,5,6,7,8,9,10)]

K_9<-cbind(K_9,sample_ID)


plot_data_K_9 <- K_9 %>% 
  gather('Population', 'Ancestry', V2:V10) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_9, aes(id, Ancestry, fill = Population)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(. ~ country, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))

  
  facet_grid(~fct_inorder(Type), switch = "x", scales = "free", space = "free")+
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))

#------------------------------------------------------------------------------------------------------------------------------#

K_10=read.table("greedy_CLUMPP/Q0.10.Q")

K_10<-cbind(K_10,sample_ID)

K_10_renamed<-left_join(K_10,id,by="id")

tbl_df(K_10_renamed)

plot_data_K_10 <- K_10_renamed %>% 
  gather('Population', 'Ancestry', V1:V10) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = Population[which.max(Ancestry)],
         assingment_prob = max(Ancestry)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

ggplot(plot_data_K_10, aes(id, Ancestry, fill = Population)) +
  geom_col() +
  facet_grid(. ~ Type, drop=TRUE, space="free", scales="free") +
  theme(panel.grid=element_blank()) +
  theme(panel.background=element_rect(fill=NA, colour="grey25")) +
  theme(panel.margin.x=grid:::unit(0.5, "lines")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(strip.text=element_text(size=12)) +
  theme(legend.key=element_rect(colour="grey25")) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
  
  
  
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", size = 4, angle = 90, hjust = .5, vjust = .5, face = "plain"))
