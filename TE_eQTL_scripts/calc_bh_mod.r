if ( length(commandArgs(TRUE)) < 1 && !exists('PREFIX')) {
  stop('Please provide a prefix.')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

if ( length(commandArgs(TRUE)) < 2 && !exists('MAXBH')) {
  MAXBH <<- 0.05
  cat('Using default maximum BH:', MAXBH, '\n')
} else if ( !exists('MAXBH') ) {
  MAXBH <<- as.numeric(commandArgs(TRUE)[2])
  cat('Using maximum BH:', MAXBH, '\n')
}

require('bigmemory')
require('Rcpp')
options(stringsAsFactors=FALSE)

data.file <- sprintf('%s_pvals.desc', PREFIX)
output.file <- sprintf('%s_bh.RData', PREFIX)

if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

Rcpp::sourceCpp('bh_mod.cpp')
bh.points <- BigMatBH(pvals@address, MAXBH)

# Returns a function that adjusts p-values according to pre-calculated bh.points
padjust.fun <- function(bh.points) {
  bh.points <- rbind(bh.points, c(NA, 1, NA))
  return( function(p) approx(bh.points$pval, bh.points$p.bh, p, method='constant', f=1)$y )
}

# Returns a function that extracts the benjamini-hochberg cut off for a given threshold
pthresh.fun <- function(bh.points) {
  return ( function (thresh) max(bh.points$pval[bh.points$p.bh <= thresh]) )
}

bh.pvals <- padjust.fun(bh.points)
bh.threshold <- pthresh.fun(bh.points)

save(bh.points, MAXBH, bh.pvals, bh.threshold, file=output.file)
