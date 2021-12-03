# LD calculation of the diversity panel 

/home/dilan/softwares/PopLDdecay/bin/PopLDdecay -InVCF /mnt/HDD-01/Dilan/SNpcallingV2/snps_fixed_GWAS_maxmiss0.5_minDP5_maf0.05_renamed.vcf.gz  -MaxDist 300  -OutStat snps_fixed_GWAS_maxmiss0.5_minDP5_maf0.05 

# Plotting 
perl  /home/dilan/softwares/PopLDdecay/bin/Plot_OnePop.pl -inFile snps_fixed_GWAS_maxmiss0.5_minDP5_maf0.05.stat.gz  -output LD_DECAY_snps_fixed_GWAS_maxmiss0.5_minDP5_maf0.05

# LD calculation for A and B genome 

# deviding vcf file in to sub genomes 
# A genome

vcftools --gzvcf quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf.gz --chr 01 --chr 02 --chr 03 --chr 04 --chr 05 --chr 06 --chr 07 --chr 08 --chr 09 --recode --recode-INFO-all --out A_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05

# B genome 
vcftools --gzvcf quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf.gz  --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --recode --recode-INFO-all --out B_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf0.05

# LD decay in sub genome A

 /home/dilan/softwares/PopLDdecay/bin/PopLDdecay -InVCF A_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05.recode.vcf -MaxDist 300  -OutStat A_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf_0.05
# LD decay in sub genome B

 /home/dilan/softwares/PopLDdecay/bin/PopLDdecay -InVCF B_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf0.05.recode.vcf  -MaxDist 300  -OutStat B_genome_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.8_minmeanDP5.maf0.05

# LD decay in different chromosomes
# seperate vcf file in to chromosomes 

for K in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 ;do vcftools --gzvcf /mnt/HDD-01/Dilan/SNpcallingV2/snps_fixed_GWAS_maxmiss0.5_minDP5_maf0.05_renamed.vcf.gz --chr $K --recode --recode-INFO-all --out ./vcf_files/$K"_chr_maxmiss0.5_minDP5_maf0.05" ; done

# LD calculation of all 
for K in ./vcf_files/*.vcf.gz; do /home/dilan/softwares/PopLDdecay/bin/PopLDdecay -InVCF $K  -MaxDist 300  -OutStat $K"_LD_decay" ; done 

# plotting 
for K in *LD_decay.stat.gz; do perl /home/dilan/softwares/PopLDdecay/bin/Plot_OnePop.pl -inFile $K --output "LD_DECAY_plots_"$K ; done 

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# selecting the genotypes only from quinoa / removing wild Chenopodium species 

vcftools --gzvcf /mnt/HDD-01/Dilan/SNpcallingV2/snps_fixed_GWAS_maxmiss0.5_minDP5_maf0.05_renamed.vcf.gz --remove wild_chenopodium.txt --recode --recode-INFO-all --out ./quinoa_diversity_set/snps_fixed_QUINOA_only_maxmiss0.5_minDP5_maf0.05

# deviding vcf file in to sub genomes 
# re-A genome
vcftools --gzvcf ./quinoa_diversity_set/snps_fixed_QUINOA_only_maxmiss0.5_minDP5_maf0.05.recode.vcf --chr 01 --chr 02 --chr 03 --chr 04 --chr 05 --chr 06 --chr 07 --chr 08 --chr 09 --recode --recode-INFO-all --out ./quinoa_diversity_set/A_genome_QUINOA_only_maxmiss0.5_minDP5_maf0.05

# re-B genome 
vcftools --gzvcf ./quinoa_diversity_set/snps_fixed_QUINOA_only_maxmiss0.5_minDP5_maf0.05.recode.vcf --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --recode --recode-INFO-all --out ./quinoa_diversity_set/B_genome_QUINOA_only_maxmiss0.5_minDP5_maf0.05



# LD decay in different chromosomes
# seperate vcf file in to chromosomes 

for K in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 ;do vcftools --gzvcf ./quinoa_diversity_set/snps_fixed_QUINOA_only_maxmiss0.5_minDP5_maf0.05.recode.vcf.gz --chr $K --recode --recode-INFO-all --out ./quinoa_diversity_set/$K"_chr_QUINOA_only_maxmiss0.5_minDP5_maf0.05" ; done

# zipping all vcf files 

gzip *.vcf

# LD calculation of all 
for K in ./quinoa_diversity_set/*.vcf.gz; do /home/dilan/softwares/PopLDdecay/bin/PopLDdecay -InVCF $K  -MaxDist 300  -OutStat $K"_LD_decay" ; done 

# plotting 
for K in *LD_decay.stat.gz; do perl /home/dilan/softwares/PopLDdecay/bin/Plot_OnePop.pl -inFile $K --output "LD_DECAY_plots_"$K ; done 



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# LD decay analysis is different echotypes

for K in highland_quinoa.txt lowland_quinoa.txt wild_chenopodium.txt; do /home/dilan/softwares/PopLDdecay/bin/PopLDdecay -InVCF /mnt/HDD-03/dilan/SNPs_V2/snpEff/snps_maxmiss0.5_minDP5_maf0.05.ANNOTATED.vcf.gz -SubPop  -MaxDist 300  -OutStat $K"_LD_decay" : done  
