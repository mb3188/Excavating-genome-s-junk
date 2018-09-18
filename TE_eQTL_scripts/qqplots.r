# This function samples k indices from 1:n starting with very dense sampling
# (every value) and then getting more and more sparse
thin <- function(n, k=2000) {
  if ( n < k ) {
    return (1:n)
  }
  
  # Figure out how many groups we need to divide n into
  # so that we end up with k values
  # We want each group to be twice as big as the previous group
  # We want to sample every index in the first group
  x <- 2
  while ( (2*x-2)/(2^x-2) > k/n ) { x <- x + 1 }
  
  # Each group is twice as big as the previous group
  end.idx <- c(round(2^(1:(x-1))*n/(2^x-2)), n)
  start.idx <- c(1, end.idx[-x]+1)
  
  # Sample the same number of values from each group
  sample.size <- ceiling(k/x)
  thin.idx <- NULL
  
  for ( i in 1:x ) {
    if ( end.idx[i] - start.idx[i] < sample.size) {
      # If there aren't enough values use all indices in group i
      thin.idx <- c(thin.idx, start.idx[i]:end.idx[i])
      # Adjust sample size accordingly
      sample.size <- ceiling((k-length(thin.idx))/(x-i))
    } else {
      # Sample from group i
      thin.idx <- c(thin.idx, sort(sample(start.idx[i]:end.idx[i], sample.size)))
    }
  }
  
  # Return the list of indices
  return ( thin.idx )
}


# Opens a new png that will contain many qq plots (80 by default)
start.image <- function(png.name, rows=8, columns=10) {
  png(png.name, width=1920, height=1440, pointsize=24)
  par(mfcol=c(rows, columns), mar=c(2, 1.8, 1, 0.2), oma=c(1, 0, 0, 0),
      mgp=c(2,0.5,0), cex=0.5, pch=19)

}

# Labels and closes png
close.image <- function(png.name) {
  mtext(png.name, cex=0.75, line=0, side=1, outer=TRUE)
  mtext(format(Sys.time(), "%Y-%m-%d %H:%M"),
        cex=0.75, line=0, side=1, adj=0.01, outer=TRUE)
  dev.off()
}


make.qq.plots <- function(columns, pvals, prefix) {

  n.genotypes <- nrow(pvals)
  thin.idx <- thin(n.genotypes)

  # Expected quantiles are the same for every phenotype
  logp.exp <- -log10(1:n.genotypes/n.genotypes)

  # Confidence intervals ala
  # Casella & Berger, 2002, 
  # 2nd edition, pg 230, Duxbury
  cint.95 <- sapply(1:n.genotypes, function (x) { qbeta(0.95, x, n.genotypes - x + 1) })
  cint.05 <- sapply(1:n.genotypes, function (x) { qbeta(0.05, x, n.genotypes - x + 1) })

  have.columns <- colnames(pvals)
  if ( is.null(have.columns) ) {
    have.columns <- 1:ncol(pvals)
  } else {
    if ( is.numeric(columns) ) {
      columns <- have.columns[columns]
    }
  }

  suffix <- paste(prefix, '.png', sep='')
  gene.num <- 0
  graph.num <- 0
  png.name <- ''

  cat('Making', length(columns), 'QQ plots..')
  for ( i in 1:length(columns) ) {
    
    gene.num <- gene.num + 1

    # Look up gene name if available
    gene.name <- columns[i]
  #   if (gene.name %in% mapped_columns ) {
  #     gene.name <- xx[which(gene.name==mapped_columns)][1]
  #     gene.name <- paste(gene.name, ' (', columns[i], ')', sep='')
  #   }

    # Make a new png file for every 80 graphs
    if ( floor((gene.num - 1)/80) + 1 != graph.num ) {
      # If graphfile counter isn't zero, we need to close previous file
      if ( graph.num != 0 ) 
        close.image(png.name)
      graph.num <- floor((gene.num - 1)/80) + 1
      
      png.name <- paste('qq', graph.num, suffix, sep='_')
      start.image(png.name)
      cat('\nNew Image\n')
    }
    
    cat(gene.num, ' ')

    # If this probe isn't in this pvalue matrix, use a placeholder
    if ( ! columns[i] %in% have.columns ) {
      plot(0, 0, type='n', bty='n',
           xlab='', ylab='',
           xaxt='n', yaxt='n',
           main=gene.name)
      mtext('No data', line=-3)
      next
    }

    # Make QQ Plot
    logp.obs <- -log10(sort(pvals[, columns[i]]))
    plot(logp.exp[thin.idx], logp.obs[thin.idx],
         pch=19, cex=0.5, 
         xlab='', ylab='', main=gene.name)

    # Add x = y diagonal line
    abline(0, 1, col='gray', lty=2)

    # Add lines for confidence intervals
    lines(logp.exp, -log10(cint.95), lty=2, col='red')
    lines(logp.exp, -log10(cint.05), lty=2, col='red')
  }

  close.image(png.name)

}