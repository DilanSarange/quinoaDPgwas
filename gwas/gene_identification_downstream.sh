# get results from the gene extraction bash script ex- genes_100kb_sorted_suggestive_DTB2.txt.bed
# get a list with gene IDs and  sli  column 7. then removed duplicates
# for example: awk -F [=-] '{print  $2}' genes_100kb_sorted_suggestive_DTF2.txt.bed | sort -u | > DTF_genes_100kb_suggestive.txt
# loop

mkdir gene_list

for K in *.bed; do awk -F [=-] '{print  $2}' $K | sort -u > ./gene_list/"list_"$K".txt" ; done 

# get gene information from the gff file
cd ./gene_list

for K in *txt; do grep -f $K /mnt/HDD-03/dilan/SNPs_V2/genes/Cq_RefSeq_v2/mRNA_quinoa_V2.gff > $K".gff"; done 

# creat bed file from the gff file

for K in *gff; do awk '{split($9,a,";");split(a[1],b,"="); print $1"\t"$4"\t"$5"\t"b[2]}' $K > $K".bed" ; done 

# creat a list of genes only with gene function and chromosome

mkdir gene_functions

for K in *gff; do awk -F "\t" '{split($9,a,";"); print $1"\t"$4"\t"$5"\t"a[1]"\t"a[8]}' $K > ./gene_functions/"functions_"$K ; done 

# creat a fasta file from the bed

mkdir fasta_files

for K in *bed; do bedtools getfasta -name $4 -fi /mnt/HDD-03/dilan/SNPs_V2/genes/Cq_RefSeq_v2/QQ74_V2_pseudomolecule.fa -bed $K -fo ./fasta_files/$K".fasta" ; done  

#extract snpEff annotations from the candidate genes 

mkdir effects_snpEff

for K in *txt; do grep -f $K /mnt/HDD-01/Dilan/final_duplicate_removed/snpEFF/maf_005/EFFECT_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minMEANDP5_0.05.ANNOTATED.vcf.gz.txt > ./effects_snpEff/"EFFECTS_snps_maxmiss0.5_minDP5_maf0.01_"$K ; done

for K in *txt; do grep -f $K /mnt/HDD-01/Dilan/final_duplicate_removed/snpEFF/maf_005/quinoa_only_snps_imputed/EFFECTS_quinoa_only_hard_filtered_snps_biallelic_FINAL_NO_duplicates_maxmiss0.5_minmeanDP5.maf_0.05.imputed_ANNOTATED.vcf.gz > ./effects_snpEff/"EFFECTS_snps_maxmiss0.5_minDP5_maf0.05_"$K ; done

# BLAST All putative candidate genes sequences
# this will take very long time therefore not recomended to run in a loop
# if nessassary undo the # and run the loop

#cd ./fasta_files
# for K in *fasta; do blastn -db nt -query $K -out $K"_BLAST_result" -perc_identity 100 -max_target_seqs 5 -outfmt 7-remote  ; done 




    