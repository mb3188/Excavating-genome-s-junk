# This script outputs an RData file with genome-wide p-values for phenotypes with BH significant hits
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

pbh.file <- sprintf('%s_pbh.RData', PREFIX)

if ( !file.exists(pbh.file) ) {
	stop('Please run bh_pvals.r script first.')
}

cat('Loading', pbh.file, '...\n')
load(pbh.file)
columns <- unique(pval.df$phenotype.idx[pval.df$p.bh < 0.05])
cat('Saving pvals for', length(columns), 'phenotypes with BH-significant hits.\n')

pvals <- pvals[, columns]

save.file <- sprintf('%s_pbh_pvals.RData', PREFIX)
save(pvals, file=save.file)