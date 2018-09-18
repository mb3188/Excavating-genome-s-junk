load.genotypes <- function(prefix) {

	# Read in genotype data from transposed PLINK format
	# plink --transposed --recode12

	gt.ids <- make.names(read.table(sprintf('%s.tfam', prefix), as.is=TRUE)[, 2])
	num.ids <- length(gt.ids)

	# SNP info in first 4 columns
	geno.data <- read.table(sprintf('%s.tped', prefix), as.is=TRUE)
	snps <- geno.data[, c(2, 1, 4)]
	colnames(snps) <- c('id', 'chrom', 'position')

	# Get rid of first 4 columns
	geno.data <- as.matrix(geno.data[, -(1:4)])

	# Extract even and odd columns (allele 1 and allele 2)
	allele.one <- geno.data[, 2*(1:num.ids)-1]
	allele.two <- geno.data[, 2*(1:num.ids)]

	# Set homozygotes to 1, heterozygotes to 0
	genotypes <- allele.one == allele.two
	# Multiply by 1 or -1 based on first allele (heterozygotes unaffected)
	genotypes <- genotypes * (2 * ((allele.one == 2) - 0.5))

	# Name dims
	colnames(genotypes) <- gt.ids
	rownames(genotypes) <- snps$id

	return ( genotypes )
}