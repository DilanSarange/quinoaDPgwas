title: "SNP based PCA analysis"
author: "Dilan Sarange"
date: "3/29/2021"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
https://github.com/DilanSarange/GWAS-workshop

## PCA analsysis using R

In this tutorial, we are going to use SNP genotypes for PCA analysis. We will use a subset of SNPs used in our quinoa GWAS analysis.  

## step 1

Open RStudio, create a new Rscript and set your working directory. 

Save your Rscript 

```{bash}
# remove "#" before running the command
# cd /mnt/HDD-08_5/GWAS_workshop/dilan
# rstudio 
```


```{r}
# working directory 
getwd()
```
## step 2

install SNPrelate package and load 

remove # before running the code in your computer 
```{r}
#if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install("SNPRelate")
library(SNPRelate)
```
## step 3

Now we will load our genotype data into r 

Both vcf and plink file format can be read using SNPrelate. But for large SNP set, plink file format is preferred because it is faster. 

Both methods are given in the following 

```{r,results='hide'}
# to use vcf file, follow these commands
# vcf.fn <- "Quinoa_SNPs_biallelic_maxmiss0.95_minmeanDP8_maf0.05_PRUNED_10k_100_r0.2.vcf.gz"
# Reformat
    #snpgdsVCF2GDS(vcf.fn, "quinoa.gds", method="biallelic.only")
    
    
# this is how we use plink genotype files    
bed.fn <-"Quinoa_SNPs_biallelic_maxmiss0.95_minmeanDP8_maf0.05_PRUNED_10k_100_r0.2.bed"
fam.fn <-"Quinoa_SNPs_biallelic_maxmiss0.95_minmeanDP8_maf0.05_PRUNED_10k_100_r0.2.fam"
bim.fn <-"Quinoa_SNPs_biallelic_maxmiss0.95_minmeanDP8_maf0.05_PRUNED_10k_100_r0.2.bim"
# this is to reformat plink file into GDS file 
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "quinoa.gds")
```

## step 4

In this step, we will run pca

```{r}
#open gds file for pca
pca_genofile<- snpgdsOpen("quinoa.gds")
# runnig the PCA 
# if you are using a large SNP data set increase the number of threads (ex: num.thread=50)
pca_quinoa_workshop <- snpgdsPCA(pca_genofile, num.thread=2)
```

## step 5

Let's extract the PCs

```{r}
# extracting first 10 PCs
quinoa_PCs <- data.frame(sample.id = pca_quinoa_workshop$sample.id,
                             PC1=pca_quinoa_workshop$eigenvect[,1], # the first eigenvector
                             PC2=pca_quinoa_workshop$eigenvect[,2], 
                             PC3=pca_quinoa_workshop$eigenvect[,3],
                             PC4=pca_quinoa_workshop$eigenvect[,4],
                             PC5=pca_quinoa_workshop$eigenvect[,5], 
                             PC6=pca_quinoa_workshop$eigenvect[,6], 
                             PC7=pca_quinoa_workshop$eigenvect[,7],
                             PC8=pca_quinoa_workshop$eigenvect[,8],
                             PC9=pca_quinoa_workshop$eigenvect[,9], 
                             PC10=pca_quinoa_workshop$eigenvect[,10],
                             stringsAsFactors = FALSE)
head(quinoa_PCs)
# write to a csv file 
write.csv(quinoa_PCs,"quinoa_PCs.csv",quote = F,row.names = F)
# variance proportion (%)
pc.percent <- pca_quinoa_workshop$varprop*100
head(round(pc.percent, 2))
```
## step 6

plotting PCs

```{r}
library(ggplot2)
ggplot(quinoa_PCs,aes(PC1,PC2))+geom_point()+
    scale_x_continuous(name="PC1 (4.52%)")+
    scale_y_continuous(name="PC2 (2.28%)")
```
© 2021 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
