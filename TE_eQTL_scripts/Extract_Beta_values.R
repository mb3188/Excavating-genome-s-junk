
#load("Expression_Methylation_bh_0.05_hits.RData")

#Exp_METH_CF_beta_values<-read.table("Exp_METH_CF_beta_values")

#BEta_values_for_25_BH_out_of_2488<-Exp_METH_CF_beta_values[,intersect(colnames(Exp_METH_CF_beta_values),TABLE_of_25_overlaping_with_24#88_and_their_interactions$phenotype)]

#write.table(BEta_values_for_25_BH_out_of_2488,"BEta_values_for_25_BH_out_of_2488")


#table_occurences_methylated_genes_NAMES_VERS<-c(t(read.table("table_occurences_methylated_genes_NAMES_VERS")))

#table_occurences_expressed_genes_NAMES_VERS<-c(t(read.table("table_occurences_expressed_genes_NAMES_VERS")))

#BEta_values_for_all_significant_interactions_1<-Exp_METH_CF_beta_values[,intersect(colnames(Exp_METH_CF_beta_values),table_occurences_methylated_genes_NAMES_VERS)]


#BEta_values_for_all_significant_interactions_FINAL<-BEta_values_for_all_significant_interactions_1[intersect(rownames(BEta_values_for_all_significant_interactions_1),table_occurences_expressed_genes_NAMES_VERS),]

#write.table(BEta_values_for_all_significant_interactions_FINAL,"BEta_values_for_all_significant_interactions_FINAL")
#save.image("BEta_ALL_significant_interactions_FINAL.RDat")


load("BEta_ALL_significant_interactions_FINAL.RDat")



library(reshape2)
reshaped_beta<-melt(as.matrix(Exp_METH_CF_beta_values))

write.table(reshaped_beta,"Reshaped_beta_ALL")


#EXTRACTED_ALL<-Exp_METH_CF_beta_values[cbind(hits.df$genotype.idx,hits.df$phenotype.idx)]

save.image("BEta_ALL_significant_interactions_FINAL.RDat")


#load("INPUT_DATA_22024_MEth_16707_GeneExp.RDat")
#pdf("HIst_METHYLATION.pdf")
#hist(as.numeric(as.vector(as.matrix(Meth_CF_mark_values_nomissing_variable))))

#dev.off()
