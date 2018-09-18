if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

if ( length(commandArgs(TRUE)) < 2 & !exists('MC.CORES')) {
  stop('Please specify number of cores on the command line or set the MC.CORES variable')
} else if ( !exists('MC.CORES') ) {
  MC.CORES <<- as.integer(commandArgs(TRUE)[2])
}

require('bigmemory')
require('biganalytics')
require('parallel')

data.file <- sprintf('%s_pvals.desc', PREFIX)

# Define dimensions of heatmap array
heatmap.probes <- 600
heatmap.snps <- 1800
heatmap <- matrix(NA, heatmap.snps, heatmap.probes)
mode(heatmap) <- 'numeric'

if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

# Split up genotypes into desired number of bins (by index)
n.genotypes <- nrow(pvals)
snp.end.idx <- round(seq(0, n.genotypes, length.out=heatmap.snps + 1)[-1])
snp.beg.idx <- c(1, snp.end.idx[-length(snp.end.idx)] + 1)

# Split up genes into desired number of bins (by index)
n.probes <- ncol(pvals)
probe.end.idx <- round(seq(0, n.probes, length.out=heatmap.probes + 1)[-1])
probe.beg.idx <- c(1, probe.end.idx[-length(probe.end.idx)] + 1)



getMean <- function(j, i) {
	rows <- snp.beg.idx[i]:snp.end.idx[i]
	cols <- probe.beg.idx[j]:probe.end.idx[j]
	return ( mean(pvals[rows, cols], na.rm=TRUE) )
}


for ( i in 1:heatmap.snps ) {
  if ( interactive() ) cat(i, 'out of', heatmap.snps, round(100*i/heatmap.snps,2), '% complete\r')
  heatmap[i, ] <- simplify2array(mclapply(1:heatmap.probes, getMean, i,
                                          mc.cores=MC.CORES))
  }

cat('\nFinished', data.file, 'Saving...\n')

# Save the heatmap data
save(heatmap, n.genotypes, n.probes, snp.end.idx, snp.beg.idx, probe.end.idx, probe.beg.idx, file=sprintf('%s_heatmap.RData', data.file))


mode(heatmap) <- 'numeric'

# Draw the heatmap image
snp.idx <- round(seq(1, n.genotypes, length.out=nrow(heatmap) + 1))
probe.idx <- round(seq(1, n.probes, length.out=ncol(heatmap) + 1))

# Approximate dimensions so that each entry in heatmap array ~ 1 pixel
png(sprintf('%s_heatmap.png', data.file), width=87+1800, height=130+600)
image(snp.idx, probe.idx, heatmap, xlab='SNPs', ylab='Genes', main=data.file)
dev.off()

