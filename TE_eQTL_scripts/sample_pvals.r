if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please provide a prefix.')
} else if ( !exists('PREFIX') ) {
  PREFIX <- commandArgs(TRUE)[1]
}

if ( length(commandArgs(TRUE)) < 2 & !exists('SAMPLE.NUM')) {
  SAMPLE.NUM <- 1e5
} else if ( !exists('PVAL.THRESHOLD') ) {
  SAMPLE.NUM <- as.numeric(commandArgs(TRUE)[2])
}


require('bigmemory')
options(stringsAsFactors=FALSE)


data.file <- sprintf('%s_pvals.desc', PREFIX)

if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

sample.idx <- sample(1.0*length(pvals), SAMPLE.NUM)

sample.row <- sample.idx %% nrow(pvals)
sample.col <- 1+(sample.idx-sample.row)/nrow(pvals)

samples <- data.frame(idx=sample.idx, row.idx=sample.row, col.idx=sample.col,
	                  pval=sapply(1:SAMPLE.NUM, function (i) pvals[sample.row[i], sample.col[i]]))


save(samples, file=sprintf('%s_samples.RData', PREFIX))
cat('Saved', nrow(samples), 'pvals.\n')
