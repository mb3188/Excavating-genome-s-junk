load("Expression_Insertion_bh_0.05_hits.RData")
load("Inserion_imput_files.RData")
phenotype_names<-phenotypes[,hits.df$genotype.idx]
hits.df$genotype.idx<-colnames(phenotype_names)


# Extract Tvmar1 drug resisince associated genes
Drug_resistance_Tvmar178_associted_B728_B7268M_genes<-hits.df[hits.df$phenotype=="Tvmar178",]

write.table(Drug_resistance_Tvmar178_associted_B728_B7268M_genes, "Drug_resistance_Tvmar178_associted_B728_B7268M_genes")


#Extract Drug resisitnce MAverick insertion associated genes 
Drug_resistance_associted_MAV88_B728_B7268M_genes<-hits.df[hits.df$phenotype=="MAV88",]

write.table(Drug_resistance_associted_MAV88_B728_B7268M_genes, "Drug_resistance_associted_MAV88_B728_B7268M_genes")

