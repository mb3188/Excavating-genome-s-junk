require('parallel')

if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please provide a prefix.')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

if ( length(commandArgs(TRUE)) < 2 & !exists('MC.CORES')) {
  MC.CORES <<- detectCores()
} else if ( !exists('MC.CORES') ) {
  MC.CORES <<- as.integer(commandArgs(TRUE)[2])
}

require('biganalytics')

data.file <- sprintf('%s_pvals.desc', PREFIX)

if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

# Calculate and store row means if not already calculated
# This can take awhile
rowavg.file <- sprintf('%s_pval_rowavgs.RData', PREFIX)
if ( file.exists(rowavg.file) ) {
  cat(rowavg.file, 'already exists. Please delete or rename if you want to rerun.\n')
} else {

  # We only do this on one core since the bottleneck is disk access
  cat('Finding average pvals for each phenotype...\n')
  pval.rowavgs <- unlist(mclapply(1:nrow(pvals), mc.cores = MC.CORES,
                         function(x) mean(pvals[x, ], na.rm = TRUE)))

  cat('Saving average pvals...\n')
  save(pval.rowavgs, file=rowavg.file)
}
