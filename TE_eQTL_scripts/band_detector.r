# This script uses the fused lasso to find 'bands'
# in rows and columns of the pvalue matrix,
# and maps the corresponding genotype/phenotype hits

if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

require('cghFLasso')
require('biganalytics')
require('bigmemory')
require('parallel')
options(stringsAsFactors=FALSE)

colavg.file <- sprintf('%s_pval_colavgs.RData', PREFIX)
if ( file.exists(colavg.file) ) {
  cat('Loading', colavg.file, '...\n')
  load(colavg.file)
} else {
  source('colavg_pvals.r')
}

flasso.col <- cghFLasso(pval.colavgs)$Esti.CopyN

# plot fused lasso estimates of column averages
png(sprintf('%s_Flasso_col.png', PREFIX), width=2000, height=600)
plot(pval.colavgs, xlab='Genes', ylab='mean pvalue',
     type='n', ylim=c(0, max(flasso.col)+0.2))
lines(flasso.col, type='s', lwd=2, col='red')
points(pval.colavgs, pch=20, col=gray(0.5), cex=1)

# Plot bonferroni hits if file exists
pbon.hits <- sprintf('%s_pbon_hits.RData', PREFIX)
if ( file.exists(pbon.hits) ) {
  cat('Loading', pbon.hits, '...\n')
  load(pbon.hits)

  # find phenotype indices of significant bonferroni hits
  phenotype <- unique(as.vector(hits.df$phenotype))
  phenotype.idx <- sapply(phenotype, function(x) unique(hits.df$phenotype.idx[which(hits.df$phenotype == x)]) )

  points(phenotype.idx, flasso.col[phenotype.idx] + 0.01,
         col='blue', bg='blue', pch=25, cex=1.5)
  text(phenotype.idx, flasso.col[phenotype.idx] + 0.02, labels=phenotype,
    cex=1.1, srt=90, adj=0, xpd=NA)
} else {
  cat('Did not find bonferroni hits file.\n')
}

dev.off()

# Calculate and store row average
# This can take 5-10 mins
#row.average <- unlist(mclapply(1: nrow(pvals), mc.cores = MC.CORES,
 #             function(x) mean(pvals[x, ], na.rm = TRUE)))
#flasso.row <- cghFLasso(row.average)$Esti.CopyN

# plot fused lasso estimates of row averages
#png(sprintf('%s_Flasso_row.png', PREFIX), width=1000, height=600)
#plot(row.average, xlab = 'SNPs', pch = 20, col = gray(0.5), ylab = 'mean pvalue', cex=1, ylim = range(flasso.row)+c(-0.05,0.15))
#abline(h = 0, col = gray(0.3))
#points(flasso.row , col = 2, type = 's', lwd = 1)
#points(genotype.idx, flasso.row[genotype.idx ], col='blue', pch = 20)
#text(genotype.idx, flasso.col[genotype.idx] + 0.1 ,labels = genotype)
#dev.off()


