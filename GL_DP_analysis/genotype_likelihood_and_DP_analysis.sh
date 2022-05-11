# extracting GT:DP:PL info from heterozygous SNPs 

for K in {10..312} ; 
do zgrep -v "^#" /mnt/HDD-02/dilan/final_duplicate_removed/quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.recode.vcf.gz |  
awk -v K="$K" '{print $9"\t"$K}' |  
awk '{split($1,a,":"); split($2,b,":"); if(a[5]=="PL") print a[1]":"a[3]":"a[4]":"a[5]"\t"b[1]"\t"b[3]"\t"b[4]"\t"b[5]; else print a[1]":"a[3]":"a[4]":"a[7]"\t"b[1]"\t"b[3]"\t"b[4]"\t"b[7] }' |  
grep '0/1' >  $K"_Heterozygous.txt" ;  
done

# extracting GT:DP:PL info from heterozygous SNPs 
for K in {10..312} ; 
 do zgrep -v "^#" /mnt/HDD-02/dilan/final_duplicate_removed/quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.recode.vcf.gz | 
 awk -v K="$K" '{print $9"\t"$K}' |  
 awk '{split($1,a,":"); split($2,b,":"); if(a[5]=="PL") print a[1]":"a[3]":"a[4]":"a[5]"\t"b[1]"\t"b[3]"\t"b[4]"\t"b[5]; else print a[1]":"a[3]":"a[4]":"a[7]"\t"b[1]"\t"b[3]"\t"b[4]"\t"b[7] }' |  
 grep '0/0' > $K"_REF_Homozygous.txt" ;
 done
 
 # extracting GT:DP:PL info from heterozygous SNPs 
 for K in {10..312} ; 
 do zgrep -v "^#" /mnt/HDD-02/dilan/final_duplicate_removed/quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.recode.vcf.gz | 
 awk -v K="$K" '{print $9"\t"$K}' |  
 awk '{split($1,a,":"); split($2,b,":"); if(a[5]=="PL") print a[1]":"a[3]":"a[4]":"a[5]"\t"b[1]"\t"b[3]"\t"b[4]"\t"b[5]; else print a[1]":"a[3]":"a[4]":"a[7]"\t"b[1]"\t"b[3]"\t"b[4]"\t"b[7] }' |  
 grep '1/1' > $K"_ALT_Homozygous.txt" ;
 done

#----------------------------------------------------------------------------------------------------------------------------#

# mean DP of the Heterozygous 
for K in {10..312}"_Heterozygous.txt" ; do awk '{sum +=$3} END {print (sum/NR)}' $K >> ./results/Hetero_DP_mean.txt ; done 

# mean PL of the heterozygous
for K in {10..312}"_Heterozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
awk '{sum +=$1} END {print (sum/NR)}' >> ./results/Hetero_PL_1_mean.txt ; done

for K in {10..312}"_Heterozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
awk '{sum +=$2} END {print (sum/NR)}' >> ./results/Hetero_PL_2_mean.txt ; done

for K in {10..312}"_Heterozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
awk '{sum +=$3} END {print (sum/NR)}' >> ./results/Hetero_PL_3_mean.txt ; done

#----------------------------------------------------------------------------------------------------------------------------#


# mean DP of the REF Homozygous
for K in {10..312}"_REF_Homozygous.txt" ; do awk '{sum +=$3} END {print (sum/NR)}' $K >> ./results/REF_Homozygous_DP_mean.txt ; done 

# mean PL of the ALT Homozygous  
for K in {10..312}"_REF_Homozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
awk '{sum +=$1} END {print (sum/NR)}' >> ./results/REF_Homozygous_PL_1_mean.txt ; done

for K in {10..312}"_REF_Homozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
awk '{sum +=$2} END {print (sum/NR)}' >> ./results/REF_Homozygous_PL_2_mean.txt ; done

for K in {10..312}"_REF_Homozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
awk '{sum +=$3} END {print (sum/NR)}' >> ./results/REF_Homozygous_PL_3_mean.txt ; done

#----------------------------------------------------------------------------------------------------------------------------#


# mean DP of the ALT Homozygous
for K in {10..312}"_ALT_Homozygous.txt" ; do awk '{sum +=$3} END {print (sum/NR)}' $K >> ./results/ALT_Homozygous_DP_mean.txt ; done 

# mean PL of the ALT Homozygous 
for K in {10..312}"_ALT_Homozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
awk '{sum +=$1} END {print (sum/NR)}' >> ./results/ALT_Homozygous_PL_1_mean.txt ; done

for K in {10..312}"_ALT_Homozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K | 
awk '{sum +=$2} END {print (sum/NR)}' >> ./results/ALT_Homozygous_PL_2_mean.txt ; done

for K in {10..312}"_ALT_Homozygous.txt" ; do awk '{split($4,a,","); print a[1]"\t"a[2]"\t"a[3]}' $K |
 awk '{sum +=$3} END {print (sum/NR)}' >> ./results/ALT_Homozygous_PL_3_mean.txt ; done

#----------------------------------------------------------------------------------------------------------------------------#
for K in {10..312}"_Heterozygous.txt" ; do awk '{sum +=$4} END {print (sum/NR)}' $K >> ./results/Hetero_GQ_mean.txt ; done
for K in {10..312}"_ALT_Homozygous.txt" ; do awk '{sum +=$4} END {print (sum/NR)}' $K >> ./results/ALT_Homozygous_GQ_mean.txt ; done
for K in {10..312}"_REF_Homozygous.txt" ; do awk '{sum +=$4} END {print (sum/NR)}' $K >> ./results/REF_Homozygous_GQ_mean.txt ; done
