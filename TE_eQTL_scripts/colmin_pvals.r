if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please provide a prefix.')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}


require('biganalytics')

data.file <- sprintf('%s_pvals.desc', PREFIX)

if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

# Check to see if we have a lot of zero p-values: did analysis complete?
percent.zero <- 100 * sum(pvals[, 1] == 0, na.rm=TRUE) / sum(!is.na(pvals[, 1]))
if ( percent.zero > 1 ) {
  cat(PREFIX, 'seems to be unfinished, aborting. (', percent.zero, '% zeros )\n')
  stop()
}

# Calculate and store column minimums if not already calculated
# This can take awhile
colmin.file <- sprintf('%s_pval_colmins.RData', PREFIX)
if ( file.exists(colmin.file) ) {
  cat(colmin.file, 'already exists. Please delete or rename if you want to rerun.\n')
} else {

  # We only do this on one core since the bottleneck is disk access
  cat('Finding min pvals for each phenotype...\n')
  pval.colmins <- colmin(pvals, na.rm=TRUE)

  cat('Saving min pvals...\n')
  save(pval.colmins, file=colmin.file)
}
