
# Required libraries
require('lrgpr')
require('parallel')
require('bigmemory') # this library is not in the default installation

# We expect at least 2 command line arguments - alternatively variables
# can be set in R prior to 'source'ing this script

# Calling from command line:
# $  Rscript glmapply_xqtl.r PREFIX MC.CORES BATCH.SIZE
# for example:
# $  Rscript glmapply_xqtl.r xqtl_data 30 5000
# will load data from xqtl_data.RData, use 30 cores and batches of 5000 genotypes

# 1 - PREFIX : if not already loaded, load data from PREFIX.RData,
#              pvalues will be stored in PREFIX_pvals.bin and PREFIX_pvals.desc
# 2 - MC.CORES : number of cores to use
# 3 - BATCH.SIZE : (optional) number of genotypes to analyze (in parallel)

if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

if ( length(commandArgs(TRUE)) < 2 & !exists('MC.CORES')) {
  
#  stop('Please specify number of cores on the command line or set the MC.CORES variable')
  MC.CORES <<- detectCores()
  } else if ( !exists('MC.CORES') ) {
  MC.CORES <<- as.integer(commandArgs(TRUE)[2])
}

if ( length(commandArgs(TRUE)) < 3 & !exists('BATCH.SIZE')) {
  # Use default
  BATCH.SIZE <<- 500
} else if ( !exists('BATCH.SIZE') ) {
  BATCH.SIZE <<- as.integer(commandArgs(TRUE)[3])
}




# If we are running via source(), data MAY already be loaded. Otherwise, load it!
if ( ! exists('genotypes') | ! exists('covars') | ! exists('phenotypes') ) {
  data.file <- sprintf('%s.RData', PREFIX)

  if ( ! file.exists(data.file) ) {
    cat('Data file not found:', data.file, '\n')
    stop('Couldn\'t find specified file.')
  }

  cat('Loading data...\n')
  load(data.file)
}

# A few checks (not exhaustive!) to see if the data is in the expected format
# if ( ! exists('genotypes') ) {stop('Missing matrix: genotypes\n')}
# if ( class(genotypes) != 'matrix' ) {stop('Not a matrix: genotypes\n')}
# if ( mode(genotypes) != 'numeric' & mode(genotypes) != 'integer' ) {stop('Not numeric: genotypes\n')}

if ( ! exists('covars') ) {stop('Missing matrix: covars\n')}
if ( ! is.null(covars) & mode(covars) != 'numeric' & mode(covars) != 'integer') {stop('Not numeric: covars\n')}

if ( ! exists('phenotypes') ) {stop('Missing matrix: phenotypes\n')}
if ( class(phenotypes) != 'matrix' ) {stop('Not a matrix: phenotypes\n')}
if ( mode(phenotypes) != 'numeric' ) {stop('Not numeric: phenotypes\n')}
if ( nrow(genotypes) != nrow(phenotypes) ) {stop('Number of rows in genotypes & phenotypes does not match\n')}

# get matrix dimensions
n <- nrow(genotypes)
n.genotypes <- ncol(genotypes)
n.phenotypes <- ncol(phenotypes)
pheno.names <- colnames(phenotypes)

# Set k to the number of columns in the design matrix (intercept + genotype + covars)
if ( is.null(covars) ) {
  # No covariates
  k <<- 2
} else if ( is.null(dim(covars)) & length(covars) == n ) {
  # One column/covariate
  k <<- 3
} else {
  k <<- 2 + ncol(covars)
}

cat('n =', n, '\nk =', k, '(intercept + genotype +', k-2, 'covariate(s))\n')
cat('n.genotypes =', n.genotypes, '\n')
cat('n.phenotypes =', n.phenotypes, '\n')

back.file <- sprintf('%s_pvals.bin', PREFIX)
desc.file <- sprintf('%s_pvals.desc', PREFIX)

# If we previously started this analysis, the big.matrix files will
# still be there and we will use those (files can be manually removed if required)
if ( file.exists(back.file) ) {
  pvals.bm <- attach.big.matrix(desc.file)
  if ( nrow(pvals.bm) != n.genotypes | ncol(pvals.bm) != n.phenotypes ) {
    stop('Found existing big.matrix files but the dimensions do not match!')
  }   
} else {
  pvals.bm <- filebacked.big.matrix(nrow=n.genotypes,
                                    ncol=n.phenotypes,
                                    type='double',
                                    backingfile=back.file,
                                    descriptorfile=desc.file,
                                    dimnames=list(NULL, pheno.names))
}


bm.desc <- describe(pvals.bm)
rm(pvals.bm) 
gc()

# Based on the number of genotypes and batch size we assign blocks of genotypes
# The last job will be smaller (if necessary)
cat('BATCH.SIZE =', BATCH.SIZE, '\n')
num.jobs <- ceiling(n.phenotypes / BATCH.SIZE)

if ( num.jobs == 1 ) {
  # If there is less than one batch, set indices to start and end.
  start.idxs <- 1
  end.idxs <- n.phenotypes
} else {
  start.idxs <- BATCH.SIZE * 0:(num.jobs - 1) + 1
  end.idxs <- c(BATCH.SIZE * 1:(num.jobs - 1), n.phenotypes)
}


# Calculate pvalues for a bunch of genotypes in parallel and save the resulting
# pvalues to the big.matrix
do.batch <- function(job.num) {
  start.idx <- start.idxs[job.num]
  end.idx <- end.idxs[job.num]

  pvals.sub <- sub.big.matrix(bm.desc,
                              firstCol=start.idx,
                              lastCol=end.idx)

  if ( pvals.sub[1, 1] != 0 && pvals.sub[nrow(pvals.sub), ncol(pvals.sub)] != 0 ) {
    cat('\nSkipping batch', job.num, 'of', num.jobs, ':', start.idx, 'to', end.idx, '...\n')

    rm(pvals.sub)
    gc()
    return ( TRUE )
  }

  
  cat('\nStarting batch', job.num, 'of', num.jobs, ':', start.idx, 'to', end.idx, '...\n')
  outer.start <- proc.time()[['elapsed']]
  
  # Calculate p-values
  Y <<- phenotypes[, start.idx:end.idx]
  if ( is.null(covars) ) {
    pvals <- glmApply2(Y ~ SNP,
                      features=genotypes,
                      nthreads=MC.CORES)$pValues
  } else {
    pvals <- glmApply2(Y ~ SNP + covars,
                      features=genotypes,
                      nthreads=MC.CORES)$pValues
  }
  
  outer.time <- proc.time()[['elapsed']] - outer.start
  cat('Completed in', outer.time, 'seconds.\n')

  cat('Saving pvals...\n')
  save.start <- proc.time()[['elapsed']]
  # Save p-values.
  pvals.sub[, ] <- pvals
  save.time <- proc.time()[['elapsed']] - save.start
  cat('Completed in', save.time, 'seconds.\n')

  rm(pvals.sub)
  gc()

  return ( TRUE )
}


# Do all the jobs!
cat('\nStarting analysis on', MC.CORES, 'cores...\n')
sapply(1:num.jobs, do.batch)
