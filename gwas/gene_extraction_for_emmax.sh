#!/bin/bash
#--------------------------------------------------------------------------------------------------------------#
#selecting SNPs that are overthe significant threshold
#suggestive threshold 8.98E-7
#significant threshold 1.8321367e-8 bonferrony correction

mkdir MAT_identification

for K in *.ps; do awk '{split($1,a,":"); print $1"\t"a[1]"\t"a[2]"\t"$3}' $K > ./MAT_identification/"MAT_plots_"$K ; done 

cd MAT_identification
#--------------------------------------------------------------------------------------------------------------------------------------------------#
# adding a header to the output file 
mkdir manhattan_plots 

for K in *.ps ; do echo -e "Marker\tChr\tPos\tp" | cat - $K > ./manhattan_plots/$K".qqman" ; done 

#--------------------------------------------------------------------------------------------------------------------------------------------------#

mkdir suggestive_SNPs
for K in *.ps; do awk '{if($4<8.98E-7) print $0}' $K > ./suggestive_SNPs/"suggestive_"$K ; done
#--------------------------------------------------------------------------------------------------------------------------------------------------#
mkdir significant_SNPs
for K in *.ps; do awk '{if($4<1.66877489e-8) print $0}' $K > ./significant_SNPs/"significant_"$K; done

rm MAT_plots_* 
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#convert SNP list to bed format and adding 50000 intervals
#suggestive + 50kb
cd suggestive_SNPs/
mkdir 50kb
for K in suggestive*; do awk '{ print "Cq"$2"\t"$3 - 50000 "\t" $3 + 50000 }' $K > ./50kb/$K".bed" ; done
cd ..
#==================================================================================================================================================#
#suggestive + 100kb
cd suggestive_SNPs/
mkdir 100kb
for K in suggestive*; do awk '{ print "Cq"$2"\t"$3 - 100000 "\t" $3 + 100000 }' $K > ./100kb/$K".bed" ; done
cd ..

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#significant + 50kb
cd significant_SNPs/
mkdir 50kb
for K in significant*; do awk '{ print "Cq"$2"\t"$3 - 50000 "\t" $3 + 50000 }' $K > ./50kb/$K".bed" ; done
#==================================================================================================================================================#
#significant + 100kb
mkdir 100kb
for K in significant*; do awk '{ print "Cq"$2"\t"$3 - 100000 "\t" $3 + 100000 }' $K > ./100kb/$K".bed" ; done
cd ..
#--------------------------------------------------------------------------------------------------------------------------------------------------#
# sorting the bed file
cd ./significant_SNPs/50kb/
for K in *bed; do sort -k1,1 -k2,2n $K > "sorted_"$K; done
cd ../../
#==================================================================================================================================================#
cd ./significant_SNPs/100kb/
for K in *bed; do sort -k1,1 -k2,2n $K > "sorted_"$K; done
cd ../../
#==================================================================================================================================================#
# sorting the bed file
cd ./suggestive_SNPs/50kb/
for K in *bed; do sort -k1,1 -k2,2n $K > "sorted_"$K; done
cd ../../
#==================================================================================================================================================#
cd ./suggestive_SNPs/100kb/
for K in *bed; do sort -k1,1 -k2,2n $K > "sorted_"$K; done
cd ../../

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#setting the gff file to be used in bedtools
#spliting the gff file and converting it to bed
# awk 'NR!=1{split($9,a,";"); print $1"\t"$4"\t"$5"\t" a[1]"\t" }' mRNA_quinoa.gff > genes.bed
#--------------------------------------------------------------------------------------------------------------------------------------------------#
# searching nearest genes
cd significant_SNPs/50kb/
mkdir genes
for K in sorted* ; do bedtools closest -a $K -b /mnt/HDD-03/dilan/SNPs_V2/genes/Cq_RefSeq_v2/genes_bed/genes_quinoa_renames_V2.bed > ./genes/"genes_50kb_"$K ; done
cd ../../
#==================================================================================================================================================#
cd significant_SNPs/100kb/
mkdir genes
for K in sorted* ; do bedtools closest -a $K -b /mnt/HDD-03/dilan/SNPs_V2/genes/Cq_RefSeq_v2/genes_bed/genes_quinoa_renames_V2.bed > ./genes/"genes_100kb_"$K ; done
cd ../../

# searching nearest genes
cd suggestive_SNPs/50kb/
mkdir genes
for K in sorted* ; do bedtools closest -a $K -b /mnt/HDD-03/dilan/SNPs_V2/genes/Cq_RefSeq_v2/genes_bed/genes_quinoa_renames_V2.bed > ./genes/"genes_50kb_"$K ; done
cd ../../
#==================================================================================================================================================#
cd suggestive_SNPs/100kb/
mkdir genes
for K in sorted* ; do bedtools closest -a $K -b /mnt/HDD-03/dilan/SNPs_V2/genes/Cq_RefSeq_v2/genes_bed/genes_quinoa_renames_V2.bed > ./genes/"genes_100kb_"$K ; done
cd ../../