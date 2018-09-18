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

# Calculate and store column minimums if not already calculated
# This can take awhile
colavg.file <- sprintf('%s_pval_colavgs.RData', PREFIX)

if ( file.exists(colavg.file) ) {
  cat(colavg.file, 'already exists. Please delete or rename if you want to rerun.\n')
} else {

  # We only do this on one core since the bottleneck is disk access
  cat('Finding average pvals for each phenotype...\n')
  # pval.colavgs <- colmean(pvals, na.rm=TRUE)
  pval.colavgs <- unlist(mclapply(1: ncol(pvals), mc.cores = MC.CORES,
                         function(x) mean(pvals[, x], na.rm = TRUE)))

  cat('Saving average pvals...\n')
  save(pval.colavgs, file=colavg.file)
}
