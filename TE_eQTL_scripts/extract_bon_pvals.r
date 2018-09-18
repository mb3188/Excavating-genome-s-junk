# This script outputs an RData file with genome-wide p-values for phenotypes with bonferonni significant hits
# These can then be used for manhattan plots (genome-wide or zoom-in)

if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

require('bigmemory')

data.file <- sprintf('%s_pvals.desc', PREFIX)

if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

minp.file <- sprintf('%s_pval_colmins.RData', PREFIX)

if ( !file.exists(minp.file) ) {
  stop('Please run colmin_pvals.r script first.')
}
cat('Loading', minp.file, '...\n')
load(minp.file)

# Figure out bonferroni cut off (don't include NAs)
p.alpha <- 0.05/ncol(pvals)/(nrow(pvals)-sum(is.na(pvals[, 1])))

# Select columns with bonferroni hits
columns <- which(pval.colmins < p.alpha)
cat('Saving pvals for', length(columns), 'phenotypes with bonferroni-significant hits.\n')

pvals <- pvals[, columns]

save.file <- sprintf('%s_pbon_pvals.RData', PREFIX)
save(pvals, file=save.file)