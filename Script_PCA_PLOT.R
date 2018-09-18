library("factoextra")
load("PCA_plotting.RData")
Tvmar1_Maverick_Insertion_noT1<-read.table("Tvmar1_Maverick_Insertion_14_strains")
genotypes_noT1<-t(Tvmar1_Maverick_Insertion_noT1)
PC_noT1<-prcomp(genotypes_noT1,scale = TRUE)
ls()
png("PCA_for_INSERTION_with_cos2_no_T1_final.png")
fviz_pca_ind(PC_noT1, 
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             col.ind = "cos2",
             repel = TRUE,
             pointsize = 4,labelsize = 6,
             ggtheme=theme(axis.text=element_text(size=14,face="bold"),
             axis.title=element_text(size=16, face="bold"),
             panel.background = element_blank()))
dev.off()
 png("PCA_for_INSERTION_without_cos2_no_T1_final.png")
fviz_pca_ind(PC_noT1, 
#col.ind = "cos2",
             repel = TRUE,
             pointsize = 4,labelsize = 6,
             ggtheme=theme(axis.text=element_text(size=14,face="bold"),
             axis.title=element_text(size=16, face="bold"),
             panel.background = element_blank()))
dev.off()
#history()
save.image("PCA_plotting.RData")