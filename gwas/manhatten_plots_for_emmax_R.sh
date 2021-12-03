echo "library(qqman)" > manhattan_plots.R

for K in *.phe.ps.qqman; do echo  $K'<-read.table("'$K'",sep="\t",header=T)' ; done >> manhattan_plots.R


for K in *.phe.ps.qqman; do echo -e 'png("'$K'.png",type="cairo",width=889,height=535,units="px")\n    manhattan('"$K"',chr="Chr",bp="Pos", snp="Marker", p="p",ylim=c(0,12), main = "'"$K"'",suggestiveline = -log10(8.98E-7),genomewideline = -log10(1.66877489e-8),cex = 1.5,cex.axis =1.3,col = c("black", "#CC3333"),chrlabs=c("A1","A2","A3","A4","A5","A6","A7","A8","A9","B1","B2","B3","B4","B5","B6","B7","B8","B9"))\n    dev.off()\n     \n    png("'$K'_QQ_plot.png",type="cairo")\n      qq('$K'$p,main = "'"$K"'")\n     dev.off() \n  #--------------------------------------------------------------------------------------# \n\n'   ; done >> manhattan_plots.R

echo "save.image()" >> manhattan_plots.R
