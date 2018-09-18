
load("Inserion_imput_files.RData")

source("regression.R") 

GenExp<-phenotypes
Insertion<-genotypes

Expression_Insertion_logp <- matrix(, nrow = ncol(GenExp), ncol =ncol(Insertion))
Expression_Insertion_beta <- matrix(, nrow = ncol(GenExp), ncol =ncol(Insertion))



for (i in 1:ncol(Insertion)){
  
    z <- data.frame(Insertion[,i])
        
    rownames(z)<-rownames(Insertion)

    fit_Expression_Methylation<-fit.lm(as.matrix(GenExp),z)

    Expression_Insertion_logp[,i] <- as.numeric(fit_Expression_Methylation$logp[,1])
    Expression_Insertion_beta[,i] <- as.numeric(fit_Expression_Methylation$beta[2,])

}


#Transform logp value into p-value and write these matrixes that were created by fitting the model

Expression_Insertion_pvals<-apply(Expression_Insertion_logp,2,function(x) 10^(-x))
rownames(Expression_Insertion_pvals)<-colnames(GenExp)
colnames(Expression_Insertion_pvals)<-colnames(Insertion)
write.table(Expression_Insertion_pvals,"Expression_Insertion_pvals")
write.table(Expression_Insertion_beta,"Expression_Insertion_beta")


library('bigmemory')

#Create big  matrix data data contain p-values and desriptior files that will be used to extract significant p#-values that are bh- adjusted 
 
x<- read.big.matrix("Expression_Insertion_pvals",type = "double", header =T,sep = " ",backingfile ="Expression_Insertion_pvals.bin", descriptor = "Expression_Insertion_pvals.desc", shared=TRUE, has.row.names = TRUE)


save.image("FINAL_finished_analysis.RDat")