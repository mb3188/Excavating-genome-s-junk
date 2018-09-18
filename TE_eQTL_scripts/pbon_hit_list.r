# This script outputs a list of bonferonni significant hits in .RData and .txt form

if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

require('bigmemory')

options(stringsAsFactors=FALSE)

# Load pvals if not already loaded
data.file <- sprintf('%s_pvals.desc', PREFIX)
if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

# Load results of colmins script (pval.colmins)
minp.file <- sprintf('%s_pval_colmins.RData', PREFIX)
if ( !file.exists(minp.file) ) {
  stop('Please run colmin_pvals.r script first.')
}
cat('Loading', minp.file, '...\n')
load(minp.file)

# Figure out bonferroni cut off (don't include NAs!)
p.alpha <- 0.05/ncol(pvals)/(nrow(pvals)-sum(is.na(pvals[, 1])))

# Select columns with bonferroni hits from pval.colmins
columns <- which(pval.colmins < p.alpha)
cat('Found', length(columns), 'phenotypes with bonferroni-significant hits.\n')

# If the input data is available, look for dimnames there
input.file <- sprintf('%s.RData', PREFIX)
if ( file.exists(input.file) ) {
  load(input.file)
  phenotype.names <- colnames(phenotypes)
  genotype.names <- colnames(genotypes)
} else {
  phenotype.names <- colnames(pvals)
  genotype.names <- NULL
}

# For the specified column (phenotype), extract all rows (genotypes) that have significant pvalues
get.col.hits <- function (col) {
  rows <- which(pvals[, col] < p.alpha)
  return( data.frame('genotype.idx'=rows, 'phenotype.idx'=col,
                     'pval'=pvals[rows, col], stringsAsFactors=FALSE) )
}

if ( length(columns) > 0 ) {
  # Get a table of hits from all the columns that contain significant hits
  hits.df <- do.call(rbind, lapply(columns, get.col.hits))

  # If we have phenotype names, add a column
  if ( ! is.null(phenotype.names) ) {
    hits.df$phenotype <- phenotype.names[hits.df$phenotype.idx]
  }

  # If we have genotype names, add a column
  if ( ! is.null(genotype.names) ) {
    hits.df$genotype <- genotype.names[hits.df$genotype.idx]
  }

  # Save RData and txt formats
  cat('Saving', nrow(hits.df), 'bonferroni-significant p-values.\n')
  save(hits.df, file=sprintf('%s_pbon_hits.RData', PREFIX))
  write.table(hits.df, sprintf('%s_pbon_hits.txt', PREFIX), row.names=FALSE)
}
