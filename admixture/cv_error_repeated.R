setwd("C:/Users/Dilan/Desktop/data_analysis_from_the_server/ADMIXTURE")
cv_error<-read.table("cv_error.txt",sep = "\t",header = T)

plot(cv_error$CV~cv_error$K)

library(ggplot2)



ggplot(data=cv_error,aes(K,CV))+geom_point()+geom_smooth()
+labs(title="Cross Validation Error",x="K",y="Cross-validation error")+
  scale_x_continuous(breaks=c(2:10),limits=c(2,10))+
  theme(plot.title = element_text(hjust = 0.5))