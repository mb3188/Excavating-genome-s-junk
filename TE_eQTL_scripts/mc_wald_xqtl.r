# We expect at least 2 command line arguments - alternatively variables
# can be set in R prior to 'source'ing this script

# Calling from command line:
# $  Rscript mc_wald_xqtl.r PREFIX MC.CORES BATCH.SIZE
# for example:
# $  Rscript mc_wald_xqtl.r xqtl_data 30 5000
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
  stop('Please specify number of cores on the command line or set the MC.CORES variable')
} else if ( !exists('MC.CORES') ) {
  MC.CORES <<- as.integer(commandArgs(TRUE)[2])
}

if ( length(commandArgs(TRUE)) < 3 & !exists('BATCH.SIZE')) {
  # Use default
  BATCH.SIZE <<- 5000
} else if ( !exists('BATCH.SIZE') ) {
  BATCH.SIZE <<- as.integer(commandArgs(TRUE)[3])
}


# Required libraries
require('compiler')
require('bigmemory') # this library is not in the default installation
library('parallel')

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
if ( ! exists('genotypes') ) {stop('Missing matrix: genotypes\n')}
if ( class(genotypes) != 'matrix' ) {stop('Not a matrix: genotypes\n')}
if ( mode(genotypes) != 'numeric' & mode(genotypes) != 'integer' ) {stop('Not numeric: genotypes\n')}

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



# Based on the number of genotypes and batch size we assign blocks of genotypes
# The last job will be smaller (if necessary)
cat('BATCH.SIZE =', BATCH.SIZE, '\n')
num.jobs <- ceiling(n.genotypes / BATCH.SIZE)

if ( num.jobs == 1 ) {
  # If there is less than one batch, set indices to start and end.
  start.idxs <- 1
  end.idxs <- n.genotypes
} else {
  start.idxs <- BATCH.SIZE * 0:(num.jobs - 1) + 1
  end.idxs <- c(BATCH.SIZE * 1:(num.jobs - 1), n.genotypes)
}


# *** do.xqtl function ********************************************************
# This function is called by do.batch() (below)

# Regress every phenotype on 1 genotype
# Requires globally defined variables: genotypes, covars, phenotypes, n, and k
do.xqtl <- cmpfun(function(gt) {

  # Make design matrix
  X <- cbind(1, genotypes[, gt], covars)
  
  # Calculate regression coefficients for all probes
  Sigma <- tryCatch(chol2inv(chol(crossprod(X))), error=function (e) { return (NULL) })

  # In case we encounter a singular matrix, return NAs
  if ( is.null(Sigma) ) {
    return ( rep(NA, ncol(phenotypes)) )
  }

  A <- tcrossprod(Sigma, X)
  Beta <- A %*% phenotypes

  residuals <- phenotypes - X %*% Beta
  
  # p-value for phenotype i
  wald.test <- function(i) {

    # Get residuals for phenotype i
    r <- residuals[, i]
    sig.sq <- sum(r^2) / (n - k)
    
    # We only care about genotype - 2nd column in design matrix
    se.b <- sqrt(sig.sq * Sigma[2, 2])
    t.stat <- Beta[2, i] / se.b
    
    # Return p-value (two-sided t-test)
    return ( 2*pt(abs(t.stat), df=(n - k), lower.tail=FALSE) )
  }
  
  # Get pvalue for each probe
  pvals <- unlist(lapply(1:ncol(phenotypes), wald.test))
  
  return ( pvals )
})


# *** do.batch function *******************************************************
# This function is called by the sapply statement below

# Calculate pvalues for a bunch of genotypes in parallel and save the resulting
# pvalues to the big.matrix
do.batch <- function(job.num) {
  require('parallel')
  start.idx <- start.idxs[job.num]
  end.idx <- end.idxs[job.num]

  # If the first AND last p-values are non-zero, skip
  if ( (pvals.bm[start.idx, 1] != 0) & (pvals.bm[end.idx, n.phenotypes] != 0) ) {
    cat('Skipping ', job.num, '. ', sep='')
    return ( FALSE )
  }

  cat('\nStarting batch', job.num, 'of', num.jobs, '...\n')
  outer.start <- proc.time()[['elapsed']]

  # Calculate p-values
  pvals <- simplify2array(mclapply(start.idx:end.idx, do.xqtl, mc.cores=MC.CORES))

  outer.time <- proc.time()[['elapsed']] - outer.start
  cat('Completed in', outer.time, 'seconds.\n')

  cat('Saving pvals...\n')
  outer.start <- proc.time()[['elapsed']]

  # Save p-values.
  pvals.sub <- sub.big.matrix(pvals.bm,
                              firstRow=start.idx,
                              lastRow=end.idx)
  # Since this is a file-backed big.matrix, this will write to disk!
  pvals.sub[, ] <- t(pvals)

  outer.time <- proc.time()[['elapsed']] - outer.start
  cat('Completed in', outer.time, 'seconds.\n')

  return ( TRUE )
}

# *** do.batch.parsave function ***********************************************
# This function is the same as do.batch except it forks a thread to save pvals
# And continues with analysis. We check to see if saving the previous batch
# completed before we fork another thread to save.

# Keep track of the process used to save
save.process <- NULL

# Calculate pvalues for a bunch of genotypes in parallel and save the resulting
# pvalues to the big.matrix


do.batch.parsave <- function(job.num) {
#  require('multicore')
  start.idx <- start.idxs[job.num]
  end.idx <- end.idxs[job.num]

  # If the first AND last p-values are non-zero, skip
  if ( (pvals.bm[start.idx, 1] != 0) & (pvals.bm[end.idx, n.phenotypes] != 0) ) {
    cat('Skipping ', job.num, '. ', sep='')
    return ( FALSE )
  }

  cat('\nStarting batch', job.num, 'of', num.jobs, '...\n')
  outer.start <- proc.time()[['elapsed']]

  # Calculate p-values
  pvals <- simplify2array(lapply(start.idx:end.idx, do.xqtl, mc.cores=MC.CORES-1))

  outer.time <- proc.time()[['elapsed']] - outer.start
  cat('Completed in', outer.time, 'seconds.\n')

  if ( !is.null(save.process) ) {
    cat('Checking if previous save finished... ')
    # collect() blocks until the save process ends
    save.time <- unlist(collect(save.process))
    cat('saving was completed in', save.time, 'seconds.\n')
  }

  if ( job.num == num.jobs ) {
    # If this is the last job, just save pvalues in this thread
    cat('Saving pvals...\n')
    save.start <- proc.time()[['elapsed']]
    # Save p-values.
    pvals.sub <- sub.big.matrix(pvals.bm,
                                firstRow=start.idx,
                                lastRow=end.idx)
    pvals.sub[, ] <- t(pvals)
    save.time <- proc.time()[['elapsed']] - save.start
    cat('Completed in', save.time, 'seconds.\n')
    return ( TRUE )
  }

  cat('Forking process to save pvals...\n')
  # Parallel forks a thread and returns a thread ID so we can check results.
  save.process <<- parallel({
    save.start <- proc.time()[['elapsed']]
    pvals.sub <- sub.big.matrix(pvals.bm,
                                firstRow=start.idx,
                                lastRow=end.idx)
    # Since this is a file-backed big.matrix, this will write to disk!
    pvals.sub[, ] <- t(pvals)
    save.time <- proc.time()[['elapsed']] - save.start
  })


  return ( TRUE )
}


# Do all the jobs!
cat('\nStarting analysis on', MC.CORES, 'cores...\n')
sapply(1:num.jobs, do.batch.parsave)