# Pick columns (phenotypes) at random and make a bunch of QQ plots

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

p.alpha <- 0.05/ncol(pvals)/(nrow(pvals)-sum(is.na(pvals[, 1])))
if ( is.null(colnames(pvals)) ) {
  columns <- which(pval.colmins < p.alpha)
} else {
  columns <- colnames(pvals)[pval.colmins < p.alpha]
}

source('qqplots.r')

make.qq.plots(columns, pvals, paste(PREFIX, '_pbon', sep=''))