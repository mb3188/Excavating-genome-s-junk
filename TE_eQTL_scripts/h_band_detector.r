# This script uses the fused lasso to find 'bands'
# in rows and columns of the pvalue matrix,
# and maps the corresponding genotype/phenotype hits

if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

# if ( length(commandArgs(TRUE)) < 2 & !exists('MC.CORES')) {
#   MC.CORES <<- 1
# } else if ( !exists('MC.CORES') ) {
#   MC.CORES <<- as.integer(commandArgs(TRUE)[2])
# }

require('cghFLasso')
require('biganalytics')
require('bigmemory')
# require('parallel')
options(stringsAsFactors=FALSE)

rowavg.file <- sprintf('%s_pval_rowavgs.RData', PREFIX)
if ( file.exists(rowavg.file) ) {
  cat('Loading', rowavg.file, '...\n')
  load(rowavg.file)
} else {
  source('rowavg_pvals.r')
}

flasso.row <- cghFLasso(pval.rowavgs)$Esti.CopyN

# plot fused lasso estimates of rowumn averages
png(sprintf('%s_Flasso_row.png', PREFIX), width=2000, height=600)
plot(pval.rowavgs, xlab='Genes', ylab='mean pvalue',
     type='n', ylim=c(0, max(flasso.row)+0.2))
lines(flasso.row, type='s', lwd=2, col='#E41A1C')
points(pval.rowavgs, pch=20, col=gray(0.5), cex=1)

# Plot bonferroni hits if file exists
pbon.hits <- sprintf('%s_pbon_hits.RData', PREFIX)
if ( file.exists(pbon.hits) ) {
  cat('Loading', pbon.hits, '...\n')
  load(pbon.hits)

  # find phenotype indices of significant bonferroni hits
  phenotype <- unique(as.vector(hits.df$phenotype))
  phenotype.idx <- sapply(phenotype, function(x) unique(hits.df$phenotype.idx[which(hits.df$phenotype == x)]) )

  points(phenotype.idx, flasso.row[phenotype.idx] + 0.01,
         col='#377EB8', bg='#377EB8', pch=25, cex=1.5)
  text(phenotype.idx, flasso.row[phenotype.idx] + 0.02, labels=phenotype,
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


