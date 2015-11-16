setwd("H:/plink-1.07-dos/")

library(car)
library(permute)
library(lattice)
library(vegan)
library(pheatmap)

#recode each marker as 0/1/2(missing NA) in plink
#--extract sig-marker.txt --keep pheno.txt --recodeA
dis_aa <- read.table("Discovery_case-only.raw", header=T)
rep_aa <- read.table("Replication_case-only.raw", header=T)
pheno <- read.table("Merged_pheno1-35_caseonly.txt")
pheno16 <- pheno[,1:18]
ab <- c("CCP-1 cit","Vim60-75 cit","Vim2-17 cit", "Fib36-52", "Fib573", "Fib 591", "CEP-1", "ptm12", "ptm13", "ptm15", "ptm35", "ptm36", "ptm37", "Fib alpha621-635 cit", "Fib alpha36-50 cit", "Fib beta60-74 cit")
names(pheno16) <- c("FID","IID",ab)
dis_test <- merge(pheno16, dis_aa, by="IID")
rep_test <- merge(pheno16, rep_aa, by="IID")

dis_markers <- c()
rep_markers <- c()
for (i in 7:33) dis_markers <- append(dis_markers,names(dis_aa[i]))
for (i in 7:22) rep_markers <- append(rep_markers,names(rep_aa[i]))

#Discovery:multi-variant logistic model
model1 <- reformulate(termlabels=dis_markers, response=ab[1])
model2 <- reformulate(termlabels=dis_markers, response=ab[2])
model3 <- reformulate(termlabels=dis_markers, response=ab[3])
model4 <- reformulate(termlabels=dis_markers, response=ab[4])
model5 <- reformulate(termlabels=dis_markers, response=ab[5])
model6 <- reformulate(termlabels=dis_markers, response=ab[6])
model7 <- reformulate(termlabels=dis_markers, response=ab[7])
model8 <- reformulate(termlabels=dis_markers, response=ab[8])
model9 <- reformulate(termlabels=dis_markers, response=ab[9])
model10 <- reformulate(termlabels=dis_markers, response=ab[10])
model11 <- reformulate(termlabels=dis_markers, response=ab[11])
model12 <- reformulate(termlabels=dis_markers, response=ab[12])
model13 <- reformulate(termlabels=dis_markers, response=ab[13])
model14 <- reformulate(termlabels=dis_markers, response=ab[14])
model15 <- reformulate(termlabels=dis_markers, response=ab[15])
model16 <- reformulate(termlabels=dis_markers, response=ab[16])
summary(lm(model1, data=dis_test))
#Filter 3 aa markers in discovery case-only lm: 27-3=24 markers remain
#AA_DRB1_11_32660115_V_P
#AA_DRB1_11_32660115_VL_A
#AA_DRB1_11_32660115_SPD_P
dis_markers_new <- dis_markers[c(-19,-21,-22)]
dis_pvalue <- c()
for (i in 1:16) {
model_new <- reformulate(termlabels=dis_markers_new, response=ab[i])
dis_pvalue <- rbind(dis_pvalue, summary(lm(model_new, data=dis_test))$coefficients[,4])
}
row.names(dis_pvalue) <- ab
write.table(dis_pvalue, "Discovery_case-only_multi-pvalue.txt", quote=F, row.names=T, col.names=T, sep="\t")
dis_sig <- ifelse(dis_pvalue < 0.05, 1, 0)
write.table(dis_sig, "Discovery_case-only_multi-sig.txt", quote=F, row.names=T, col.names=T, sep="\t")
#Heatmap based on multi-variant logistic regression significance jaccard distance
dis_sig_new <- dis_sig[,-1]
dis_dist <- vegdist(dis_sig_new, method="jaccard") 
#warning:missing values in jaccard distance results
pheatmap(dis_sig_new, scale="none", clustering_distance_col=dis_dist, cluster_rows=T, cluster_cols=F, fontsize=12, fontsize_rows=5)
#multi-collinearity
vif(lm(model1, data=dis_test))

#Replication:multi-variant logistic model
model1 <- reformulate(termlabels=rep_markers, response=ab[1])
model2 <- reformulate(termlabels=rep_markers, response=ab[2])
model3 <- reformulate(termlabels=rep_markers, response=ab[3])
model4 <- reformulate(termlabels=rep_markers, response=ab[4])
model5 <- reformulate(termlabels=rep_markers, response=ab[5])
model6 <- reformulate(termlabels=rep_markers, response=ab[6])
model7 <- reformulate(termlabels=rep_markers, response=ab[7])
model8 <- reformulate(termlabels=rep_markers, response=ab[8])
model9 <- reformulate(termlabels=rep_markers, response=ab[9])
model10 <- reformulate(termlabels=rep_markers, response=ab[10])
model11 <- reformulate(termlabels=rep_markers, response=ab[11])
model12 <- reformulate(termlabels=rep_markers, response=ab[12])
model13 <- reformulate(termlabels=rep_markers, response=ab[13])
model14 <- reformulate(termlabels=rep_markers, response=ab[14])
model15 <- reformulate(termlabels=rep_markers, response=ab[15])
model16 <- reformulate(termlabels=rep_markers, response=ab[16])
summary(lm(model1, data=rep_test))
#Filter 5 aa markers in replication case-only lm: 16-5=11 markers remain
#AA_DRB1_11_32660115_SD_P
#AA_DRB1_11_32660115_VL_A
#AA_DRB1_11_32660115_VD_P
#AA_DRB1_11_32660115_SPD_P
#AA_DRB1_11_32660115_SLD_P
rep_markers_new <- rep_markers[c(-12,-13,-14,-15,-16)]
rep_pvalue <- c()
for (i in 1:16) {
  model_new <- reformulate(termlabels=rep_markers_new, response=ab[i])
  rep_pvalue <- rbind(rep_pvalue, summary(lm(model_new, data=rep_test))$coefficients[,4])
}
row.names(rep_pvalue) <- ab
write.table(rep_pvalue, "Replication_case-only_multi-pvalue.txt", quote=F, row.names=T, col.names=T, sep="\t")
rep_sig <- ifelse(rep_pvalue < 0.05, 1, 0)
write.table(rep_sig, "Replication_case-only_multi-sig.txt", quote=F, row.names=T, col.names=T, sep="\t")
#Heatmap based on multi-variant logistic regression significance jaccard distance
rep_sig_new <- rep_sig[-10,-1]
rep_dist <- vegdist(rep_sig_new, method="jaccard") 
#warning:missing values in jaccard distance results
pheatmap(rep_sig_new, scale="none", clustering_distance_col=rep_dist, cluster_rows=T, cluster_cols=F, fontsize=12, fontsize_rows=5)
#multi-collinearity
vif(lm(model1, data=rep_test))
