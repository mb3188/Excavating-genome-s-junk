# Make QQ plots for all benjamini-hochberg significant hits

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

source('qqplots.r')

make.qq.plots(columns, pvals, paste(PREFIX, '_pbh', sep=''))