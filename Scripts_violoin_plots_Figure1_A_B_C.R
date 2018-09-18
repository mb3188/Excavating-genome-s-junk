########### Plot GENES WITH TE vs GENES without TE regardless of up or down or both sides

par(cex.axis=1.25, lwd=3,mfrow=c(1,1))

#png(filename="All_TE_Presence_Absence_G3_expression")
#par(cex.axis=1.25, lwd=3)
vioplot(c(Average_Genes_TE_both_sides_expression_G3, Average_TE_Upstream_of_GENE_FINAL_list_expression_G3, Average_TE_Downstream_of_GENE_FINAL_list_expression_G3), Average_GENES_without_TE_insertion_FINAL_list_expression_G3, names=c("TE", "NO TE"))
vioplot(c(Average_Genes_TE_both_sides_expression_G3, Average_TE_Upstream_of_GENE_FINAL_list_expression_G3, Average_TE_Downstream_of_GENE_FINAL_list_expression_G3),at=1,col="red",add=TRUE)
vioplot(Average_GENES_without_TE_insertion_FINAL_list_expression_G3,at=2,col="blue",add=TRUE)
axis(side = 2, font = 2)
axis(side = 1, font = 1)
#dev.off()


## TEST DISTANCE VS EXPRESSION 
require(stats)
reg<-lm(abs(Combined_expression_all_genes_with_TE_distance$V2)~rowMeans(Combined_expression_all_genes_with_TE_insertion))

reg<-lm(rowMeans(Combined_expression_all_genes_with_TE_insertionabs)~abs(Combined_expression_all_genes_with_TE_distance$V2))


coeff=coefficients(reg)
# equation of the line : 
eq = paste0("y = ", round(coeff[2],1), "*x ", round(coeff[1],1))
# plot
png("Distance_vs_Expression_all_TE_genes_G3")
par(cex.axis=1.25, lwd=3)
plot(abs(Combined_expression_all_genes_with_TE_distance$V2)~rowMeans(Combined_expression_all_genes_with_TE_insertion), main=eq)
abline(reg, col="red") 
axis(side = 2, font = 2)
axis(side = 1, font = 1)
dev.off()


####test for difference between BspAs

#Wilcoxon_Test_BspA<-wilcox.test(BspA_overlapp_with_TE_list_expression_table_G3_expression, Expressed_genes_NO_TE_BspA_list_expression_table_G3_expression, alternative = "two.sided")

#png("BspA_TE_noTE_expression_G3.png")
#par(cex.axis=1.25, lwd=3)
vioplot(BspA_overlapp_with_TE_list_expression_table_G3_expression, Expressed_genes_NO_TE_BspA_list_expression_table_G3_expression, names=c("TE", "NO TE"))
vioplot(BspA_overlapp_with_TE_list_expression_table_G3_expression,at=1,col="red",add=TRUE)
vioplot(Expressed_genes_NO_TE_BspA_list_expression_table_G3_expression,at=2,col="blue",add=TRUE)
axis(side = 2, font = 2)
axis(side = 1, font = 1)
#dev.off()

