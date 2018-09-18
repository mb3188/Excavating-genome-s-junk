if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please provide a prefix.')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

if ( length(commandArgs(TRUE)) < 2 & !exists('MC.CORES')) {
  MC.CORES <<- 1
} else if ( !exists('MC.CORES') ) {
  MC.CORES <<- as.integer(commandArgs(TRUE)[2])
}

if ( length(commandArgs(TRUE)) < 3 & !exists('PVAL.THRESHOLD')) {
  PVAL.THRESHOLD <<- NULL
} else if ( !exists('PVAL.THRESHOLD') ) {
  PVAL.THRESHOLD <<- as.numeric(commandArgs(TRUE)[3])
}


require('bigmemory')
require('biganalytics')
require('parallel')
options(stringsAsFactors=FALSE)


# Function to find raw p-value cutoffs for BH thresholds
find.pbh.boundary <- function(pvals, p.bh, threshold=0.05) {
  if ( length(threshold) > 1 ) {
    return ( sapply(threshold, function (t) find.boundary(pvals, p.bh, t)) )
  }
  stopifnot(length(pvals)==length(p.bh))
  
  if ( min(p.bh) > threshold || max(p.bh) < threshold ) {
    return ( NA )
  }
  lower.bound <- max(pvals[p.bh <= threshold])
  upper.bound <- min(pvals[p.bh >= threshold])
  boundary <- mean(lower.bound, upper.bound)
  if ( boundary > upper.bound || boundary < lower.bound ) {
    return ( NA )
  }
  num.digits <- ceiling(-log10(upper.bound-lower.bound))
  approx.boundary <- round(boundary, num.digits)

  # Reduce number of digits until we are no longer in bounds
  while ( approx.boundary <= upper.bound && approx.boundary > lower.bound ) {
    num.digits <- num.digits - 1
    approx.boundary <- round(boundary, num.digits)
  }
  # Return rounded value
  return ( round(boundary, num.digits+1) )
}





data.file <- sprintf('%s_pvals.desc', PREFIX)

if ( !exists('pvals') ) {
  cat('Loading', data.file, '...\n')
  pvals <- attach.big.matrix(data.file)
}


# Calculate and store column minimums if not already calculated
# This can take awhile
colmin.file <- sprintf('%s_pval_colmins.RData', PREFIX)
if ( file.exists(colmin.file) ) {
  cat('Loading min pvals...\n')
  load(colmin.file)
} else {
  cat('Please run script colmin_pvals.r to find column minimums before running this script...\n')

  # source('colmin_pvals.r')
}

# If threshold undefined (NULL), estimate threshold so that we select ~ 2.5% of columns
if ( is.null(PVAL.THRESHOLD) ) {
  PVAL.THRESHOLD <- signif(quantile(pval.colmins, 0.025), 2)
} 

cat('Using pvalue threshold', PVAL.THRESHOLD, 'and', MC.CORES, 'cores.\n')

# Make dataframe to store smallest pvals
pval.df <- data.frame()
phenotype.idx <- 1:ncol(pvals)
cat('Selecting columns...\n')
use.columns <- which(pval.colmins < PVAL.THRESHOLD)

cat('Pulling pvals from', length(use.columns), 'columns...\n')
get.pvals <- function(i) {
  idx <- mwhich(pvals, i, PVAL.THRESHOLD, comps='lt')
  return ( data.frame('phenotype.idx'=phenotype.idx[i],
                      'genotype.idx'=idx,
                      'pval'=pvals[idx, i]) )
}

pval.df <- do.call(rbind, mclapply(use.columns, get.pvals, mc.cores=MC.CORES))

cat('Calculating adjusted pvals...\n')
p.bh <- p.adjust(as.numeric(pval.df$pval), method='BH', n=length(pvals))
pval.df$p.bh <- p.bh

cat('Adjusted pval range:', range(p.bh), '\n')

save(pval.df, file=sprintf('%s_pbh.RData', PREFIX))
cat('Saved', length(p.bh), 'pvals.\n')

if ( max(p.bh) < 0.05 ) {
  cat('\nWARNING:\n')
  cat('The largest extracted p-value is Benjamini-Hochberg significant, please rerun script using a larger p-value cutoff, e.g.\n')
  cat('        Rscript bh_pvals.r', MC.CORES, PVAL.THRESHOLD*10, '\n')

  stop('Incomplete list of hits with current threshold.\n')
} else if ( min(p.bh) > 0.05 ) {
  stop('No Benjamini-Hochberg significant hits detected!')
} else {
  cat('Threshold for Benjamini-Hochberg significance (0.05) is:',
       find.pbh.boundary(pval.df$pval, p.bh, 0.05), '\n')
}

# Count phenotypes with Benjamini-Hochberg hits
num.hits <- unique(pval.df$phenotype.idx[pval.df$p.bh < 0.05])
cat('Found', length(num.hits), 'phenotypes with Benjamini-Hochberg-significant hits.\n')

if ( nrow(pval.df) > 0 ) {
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

  # Get the significant hits from pval.df
  hits.df <- pval.df[p.bh < 0.05, ]

  # Rename columns to match pbon_hit_list.r output
  colnames(hits.df) <- c('phenotype.idx', 'genotype.idx', 'pval', 'p.bh')

  # If we have phenotype names, add a column
  if ( ! is.null(phenotype.names) ) {
    hits.df$phenotype <- phenotype.names[hits.df$phenotype.idx]
  }

  # If we have genotype names, add a column
  if ( ! is.null(genotype.names) ) {
    hits.df$genotype <- genotype.names[hits.df$genotype.idx]
  }

  # Save RData and txt formats
  cat('Saving', nrow(hits.df), 'Benjamini-Hochberg-significant hits.\n')
  save(hits.df, file=sprintf('%s_pbh_hits.RData', PREFIX))
  write.table(hits.df, sprintf('%s_pbh_hits.txt', PREFIX), row.names=FALSE)
}
