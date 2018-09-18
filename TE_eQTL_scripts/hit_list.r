# This script outputs a list of significant hits in .RData and .txt form


require('bigmemory')
options(stringsAsFactors=FALSE)


if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

if ( length(commandArgs(TRUE)) < 2 & !exists('UBTYPE')) {
  UBTYPE <<- 'bon'
  cat('Using default upper bound type:', UBTYPE, '\n')
} else if ( !exists('UBTYPE') ) {
  UBTYPE <<- commandArgs(TRUE)[2]
  if ( ! UBTYPE %in% c('bon', 'bh', 'raw') )
    stop(paste('Upper bound type "', UBTYPE, '" not recognized. Please specify a upper bound type of bon, bh, or raw.', sep=''))
  cat('Using upper bound type:', UBTYPE, '\n')
}

if ( length(commandArgs(TRUE)) < 3 & !exists('UBOUND')) {
  UBOUND <<- 0.05
  cat('Using default upper bound:', UBOUND, '\n')
} else if ( !exists('UBOUND') ) {
  UBOUND <<- as.numeric(commandArgs(TRUE)[3])
  cat('Using upper bound:', UBOUND, '\n')
}

if ( length(commandArgs(TRUE)) < 4 & !exists('LBTYPE')) {
  LBTYPE <<- NULL
  LBOUND <<- NULL
} else if ( !exists('LBTYPE') ) {
  LBTYPE <<- commandArgs(TRUE)[4]
  if ( ! UBTYPE %in% c('bon', 'bh', 'raw') )
    stop(paste('Lower bound type "', LBTYPE, '" not recognized. Please specify a lower bound type of bon, bh, or raw.', sep=''))
  cat('Using lower bound type:', LBTYPE, '\n')

  if ( length(commandArgs(TRUE)) < 5 & !exists('LBOUND')) {
    LBOUND <<- 0.05
    cat('Using default lower bound:', LBOUND, '\n')
  } else if ( !exists('LBOUND') ) {
    LBOUND <<- as.numeric(commandArgs(TRUE)[5])
    cat('Using lower bound:', LBOUND, '\n')
  }
}

output.name <- paste(PREFIX, UBTYPE, UBOUND, sep='_')
if ( !is.null(LBTYPE) ) {
  output.name <- paste(output.name, 'to', LBTYPE, LBOUND, sep='_')
}

# Load pvals if not already loaded
data.file <- sprintf('%s_pvals.desc', PREFIX)
if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}

# Load results of colmins script (pval.colmins)
minp.file <- sprintf('%s_pval_colmins.RData', PREFIX)
if ( !file.exists(minp.file) ) {
  cat('Calculating column minimums...\n')
  source('colmin_pvals.r')
} else {
  cat('Loading', minp.file, '...\n')
  load(minp.file)
}

# Load BH cutoffs if needed
if ( UBTYPE == 'bh' || (!is.null(LBTYPE) && LBTYPE == 'bh') ) {
  bh.file <- sprintf('%s_bh.RData', PREFIX)
  if ( !file.exists(bh.file) ) {
  cat('Calculating BH cutoffs...\n')
  MAXBH <- UBOUND
  source('calc_bh.r')
} else {
  cat('Loading', bh.file, '...\n')
  load(bh.file)
  }

}

if ( UBTYPE == 'bon' ) {
  p.upper <- UBOUND/length(pvals)
} else if ( UBTYPE == 'bh' ) {
  p.upper <- max(bh.points$pval[bh.points$p.bh <= UBOUND])
} else {
  p.upper <- UBOUND
}

if ( is.null(LBTYPE) ) {
  p.lower <- 0
} else if ( LBTYPE == 'bon' ) {
  p.lower <- LBOUND/length(pvals)
} else if ( LBTYPE == 'bh' ) {
  p.lower <- max(bh.points$pval[bh.points$p.bh <= LBOUND])
} else {
  p.lower <- LBOUND
}

cat('\nExtracting pvals with (', p.lower, ' <= p < ', p.upper, ')\n\n', sep='')


# Select columns with possible hits from pval.colmins
columns <- which(pval.colmins < p.upper)
cat('Scanning', length(columns), 'phenotypes with possible hits.\n')

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
  rows <- which(pvals[, col] < p.upper & pvals[, col] >= p.lower)
  if ( length(rows) == 0 ) return ( NULL )
  return( data.frame('genotype.idx'=rows, 'phenotype.idx'=col,
                     'pval'=pvals[rows, col], stringsAsFactors=FALSE) )
}

if ( length(columns) > 0 ) {
  
  require('parallel')
  # Get a table of hits from all the columns that contain significant hits
  hits.df <- do.call(rbind, mclapply(columns, get.col.hits, mc.cores=detectCores()))

  cat('Found', nrow(hits.df), 'pvalues to save.\n')

  hits.df <- hits.df[order(hits.df$pval), ]

  # If we have phenotype names, add a column
  if ( ! is.null(phenotype.names) ) {
    hits.df$phenotype <- phenotype.names[hits.df$phenotype.idx]
  }

  # If we have genotype names, add a column
  if ( ! is.null(genotype.names) ) {
    hits.df$genotype <- genotype.names[hits.df$genotype.idx]
  }
  
  # If we have BH points calculated, add p.bh column
  if ( exists('bh.pvals') ) {
    hits.df$p.bh <- bh.pvals(hits.df$pval)
  }

  # Save RData and txt formats
  cat('Saving ', nrow(hits.df), ' p-values (', 100*nrow(hits.df)/length(pvals),' percent of total values).\n', sep='')

  attr(hits.df, 'upper.bound') <- p.upper
  attr(hits.df, 'upper.bound.type') <- paste(UBTYPE, UBOUND)
  attr(hits.df, 'lower.bound') <- p.lower
  attr(hits.df, 'lower.bound.type') <- paste(LBTYPE, LBOUND)

  save(hits.df, file=sprintf('%s_hits.RData', output.name))
  write.table(hits.df, sprintf('%s_hits.txt', output.name), row.names=FALSE)
}
