#' Fit linear model
#'
#' @param expr Numeric matrix of expression values with samples in rows and
#' genes in columns. Row names must match row names of sample.data data.frame.
#' @param sample.data Data.frame containing covariates which should be
#' jointly fit in the model. Will be encoded using model.matrix function. Row
#' names must match row names of expr matrix.
#' @param save.input Whether to include input as part of result object
#'
#' @return A list object with results of model fit, including
#' \describe{
#'   \item{\code{logp}}{negative log base 10 p values for each gene for each covariate in sample.data}
#'   \item{\code{beta}}{regression coefficients for each gene for each encoded covariate}
#'   \item{\code{y}}{gene expression matrix for samples in model}
#'   \item{\code{z}}{design/model matrix with encoded covariates}
#'   \item{\code{df.map}}{Indices of z/beta matrix columns/rows which correspond to a
#'   single column in sample.data (otherwise NA)}
#'   \item{\code{z.map}}{Indices of sample.data/logp columns which correspond to each
#'   column of z}
#'   \item{\code{input}}{(Optional) if specified, the input is saved here.}
#' }
#'
#' @importFrom stats model.matrix pf pt
#' @export
fit.lm <- function(expr, sample.data, save.input=FALSE) {
  stopifnot(rownames(expr)==rownames(sample.data))

  # Generate model matrix
  z <- model.matrix(~., data=sample.data)

  k <- ncol(z)
  n <- nrow(z)

  if ( n < nrow(expr) ) {
    message(nrow(expr)-n, ' samples dropped due to missing data.')
  }

  y <- expr[rownames(z), , drop=FALSE]
  # Solve for sigma matrix
  sigma <- chol2inv(chol(crossprod(z)))

  # Intermittent issue causes NAs in inverted matrix, usually rerunning fixes it.
  stopifnot(sum(is.na(sigma)) == 0)

  # Calculate "hat" matrix (a) & regression coefficients (beta)
  a <- tcrossprod(sigma, z)
  beta <- a %*% y

  # Calculate t-statistics for all variables in model
  resid <- y - z %*% beta
  sig.sq <- apply(resid^2, 2, sum)/(n-k)
  se.b <- do.call(c, lapply(1:k, function (i) sqrt(sig.sq * sigma[i, i])))
  dim(se.b) <- c(ncol(y), k)
  t.stat <- t(beta) / se.b

  # Calculate log p for t-statistics
  zlogp <- -pt(abs(t.stat), df=(n-k), lower.tail=FALSE, log.p=TRUE)/log(10)-log10(2)
  colnames(zlogp) <- colnames(z)
  rownames(zlogp) <- colnames(y)

  # Map data columns to model columns
  z.map <- attr(z, 'assign')
  tzm <- table(z.map[-1])
  single.df.idx <- as.numeric(names(tzm)[tzm==1])
  single.z.idx <- which(z.map %in% single.df.idx)

  # Reorganize log p values by data columns
  logp <- matrix(NA, nrow=ncol(y), ncol=ncol(sample.data))
  dimnames(logp) <- list(colnames(expr), colnames(sample.data))
  logp[, single.df.idx] <- zlogp[, single.z.idx]

  multi.df.idx <- which(! 1:ncol(sample.data) %in% single.df.idx)

  # Use F test to fill in log p values for categorical variables
  for ( i in multi.df.idx ) {
    logp[, i] <- lm.ftest(y, z, z[, z.map != i], resid)
  }

  rownames(beta) <- colnames(z)
  df.map <- rep(NA, ncol(sample.data))
  df.map[single.df.idx] <- single.z.idx

  names(z.map) <- colnames(z)
  names(df.map) <- colnames(sample.data)

  fit <- list()
  fit$logp <- logp
  fit$beta <- beta
  fit$y <- y
  fit$z <- z
  fit$df.map <- df.map
  fit$z.map <- z.map
  if ( save.input ) {
    fit$input <- list()
    fit$input$expr <- expr
    fit$input$sample.data <- sample.data
  }

  return (bh.adjust(fit))
}


# Calculate F test pvalue for multiple response variables
lm.ftest <- function(y, z, z0, resid.full) {
  if ( is.null(dim(y)) )
    y <- matrix(y, ncol = 1)

  n <- nrow(y)

  if ( n != nrow(z) ) {
    stop('rows of z must match rows of y.')
  }

  if ( missing(z0) ) {
    z0 <- matrix(1, nrow = n, ncol=n)
  }

  k <- ncol(z)
  k0 <- ncol(z0)

  # Don't recalculate residuals for full model if we already have them
  if ( missing(resid.full) ) {
    sigma <- chol2inv(chol(crossprod(z)))
    if ( is.null(sigma) ) {
      return ( rep(0, ncol(y)) )
    }
    a <- tcrossprod(sigma, z)
    beta <- a %*% y
    resid.full <- y - z %*% beta
  }

  # Calculate residuals for null model
  sigma0 <- chol2inv(chol(crossprod(z0)))
  a0 <- tcrossprod(sigma0, z0)
  beta0 <- a0 %*% y
  resid.null <- y - z0 %*% beta0

  # Calculate F statistic
  sse.full<- apply(resid.full^2, 2, sum)
  sse.null <- apply(resid.null^2, 2, sum)
  msr <- (sse.null - sse.full)/(k-k0)
  mse <- sse.full/(n-k)
  f.stat <- msr/mse

  # Return log10 scaled pvalue for each gene
  logp <- -pf(f.stat, df1=k-k0, df2=n-k, lower.tail=FALSE, log.p=TRUE)/log(10)

  return ( logp )
}


# Calculate likelihood ratio test pvalue for multiple response variables
lm.lrtest <- function(y, z, z0, resid.full) {
  if ( is.null(dim(y)) )
    y <- matrix(y, ncol = 1)
  
  n <- nrow(y)
  
  if ( n != nrow(z) ) {
    stop('rows of z must match rows of y.')
  }
  
  if ( missing(z0) ) {
    z0 <- matrix(1, nrow = n, ncol=n)
  }
  
  k <- ncol(z)
  k0 <- ncol(z0)
  
  # Don't recalculate residuals for full model if we already have them
  if ( missing(resid.full) ) {
    sigma <- chol2inv(chol(crossprod(z)))
    if ( is.null(sigma) ) {
      return ( rep(0, ncol(y)) )
    }
    a <- tcrossprod(sigma, z)
    beta <- a %*% y
    resid.full <- y - z %*% beta
  }
  
  # Calculate residuals for null model
  sigma0 <- chol2inv(chol(crossprod(z0)))
  a0 <- tcrossprod(sigma0, z0)
  beta0 <- a0 %*% y
  resid.null <- y - z0 %*% beta0
  
  # Calculate F statistic
  sse.full<- apply(resid.full^2, 2, sum)
  sse.null <- apply(resid.null^2, 2, sum)
  var.full <- apply(resid.full, 2, var)
  var.null <- apply(resid.null, 2, var)
  
  loglik.full <- -(n/2)*log(2*pi)-(n/2)*log(var.full)-sse.full/(2*var.full)
  loglik.null <- -(n/2)*log(2*pi)-(n/2)*log(var.null)-sse.null/(2*var.null)
  
  lr.stat <- 2*(loglik.full-loglik.null)
  
  # Return log10 scaled pvalue for each gene
  logp <- -pchisq(lr.stat, df=k-k0, lower.tail=FALSE, log.p=TRUE)/log(10)
  
  return ( logp )
}


# Add Benjamini-Hochberg correction information to fitted model
bh.adjust <- function(fit, alpha=0.05) { 
  fit$bh.threshold <- apply(fit$logp, 2, bh.threshold, alpha)
  return ( fit )
}

bh.threshold <- function(logp, alpha=0.05) {
  pvals <- 10^-logp
  bh.pvals <- p.adjust(pvals, method='BH')
  bh.logp <- -log10(bh.pvals)
  
  # Smallest significant logp
  lower.bound <- max(logp[bh.pvals>=alpha])
  # Largest non-significant logp
  upper.bound <- suppressWarnings(min(logp[bh.pvals<alpha]))
  
  if ( is.na(lower.bound) | is.na(upper.bound) )
    return ( -log10(0.05/length(pvals)) )
  
  between(lower.bound, upper.bound)
}


# Return a roundish number between the threshold values
between <- function(lower, upper) {
  digits <- 0
  floord <- function (x, d) floor(10^d*x)/10^d 
  while ( floord(upper, digits) <= lower ) 
    digits <- digits + 1
  floord(upper, digits)
}


# Convenience function for calculating additional F tests for fitted model
ftest <- function(fit, terms) {
  
  if ( all(terms %in% colnames(fit$logp)) ) {
    z.terms <- match(match(terms, colnames(fit$logp)), fit$z.map)
  } else if ( all(terms %in% colnames(fit$z)) ) {
    z.terms <- match(terms, colnames(fit$z))
  } else {
    z.terms <- terms
  }
    lm.ftest(fit$y, fit$z, fit$z[, -z.terms], fit$y-fit$z %*% fit$beta)
}

# Convenience function for calculating additional LR tests for fitted model
lrtest <- function(fit, terms) {
  
  if ( all(terms %in% colnames(fit$logp)) ) {
    z.terms <- match(match(terms, colnames(fit$logp)), fit$z.map)
  } else if ( all(terms %in% colnames(fit$z)) ) {
    z.terms <- match(terms, colnames(fit$z))
  } else {
    z.terms <- terms
  }
  lm.lrtest(fit$y, fit$z, fit$z[, -z.terms], fit$y-fit$z %*% fit$beta)
}


# Estimate best "direction" for composite tests
comp.sign <- function(fit, terms) {
  if ( all(terms %in% colnames(fit$logp)) ) {
    p.terms <- match(terms, colnames(fit$logp))
    z.terms <- match(p.terms, fit$z.map)
  } else if ( all(terms %in% colnames(fit$z)) ) {
    z.terms <- match(terms, colnames(fit$z))
    p.terms <- match(z.terms, fit$df.map)
  }
  
  # Initialize with signs of first term
  stopifnot(length(z.terms) > 0)
  csign <- sign(fit$beta[z.terms[1], ])
  clogp <- fit$logp[, p.terms[1]]
  
  # Loop through the rest of the terms updating sign as required
  if ( length(z.terms) > 1 ) {
    for ( i in 2:length(z.terms) ) {
      zt <- z.terms[i]
      pt <- p.terms[i]
      
      cor.sign <- sign(cor(fit$z[, z.terms[1]], fit$z[, zt]))
      
      # If logp is greater, use sign of this term's coefficient
      update.sign <- clogp < fit$logp[, pt]
      if ( any(update.sign) ) {
        csign[update.sign] <- cor.sign * sign(fit$beta[zt, update.sign])
        clogp[update.sign] <- fit$logp[update.sign, pt]
      }
    }
  }
  return ( csign ) 
}


# Plot logp vs logp to compare models 
pvp.plot <- function(x, y, x.lab='', y.lab='',
                     x.thresh, y.thresh,
                     x.thresh.lab=expression(p[BH]==alpha), y.thresh.lab=expression(p[BH]==alpha),
                     add.legend=TRUE, filename, ...) {
  if ( missing(x.thresh) ) 
    x.thresh <- bh.threshold(abs(x))
  if ( missing(y.thresh) )
    y.thresh <- bh.threshold(abs(y))
  
  col.code <- rep(1, length(x))
  col.code[abs(x) > x.thresh & abs(y) > y.thresh] <- 2
  col.code[abs(x) > x.thresh & abs(y) < y.thresh] <- 3
  col.code[abs(x) < x.thresh & abs(y) > y.thresh] <- 4
  col.code[((x > 0 & y < 0) | (x < 0 & y > 0)) & (abs(x) > x.thresh | abs(y) > y.thresh)] <- 5
  colors <- RColorBrewer::brewer.pal(5, 'Set1')
  
  if ( !missing(filename) ) {
    png(paste0('results/', filename), width=800, height=600, pointsize=20)
    par(pch=20, mar=c(3, 3, 2, 2), mgp=c(1.5, 0.5, 0), bty='n', las=1)
  }
  plot(x, y, xlab=x.lab, ylab=y.lab, col=colors[col.code], ...)
  
  abline(h=y.thresh, v=x.thresh, lty=3)
  abline(h=bon.logp, v=bon.logp, col='#cccccc')
  axis(3, c(x.thresh, bon.logp),
       c(x.thresh.lab, expression(p[Bon]==alpha)),
       lwd=0, mgp=c(0, 0, 0), cex.axis=0.6)
  axis(4, c(y.thresh, bon.logp),
       c(y.thresh.lab, expression(p[Bon]==alpha)),
       lwd=0, mgp=c(0, 0, 0), cex.axis=0.6)
  if ( any(x < 0 | y < 0) ) {
    abline(h=-y.thresh, v=-x.thresh, lty=3)
    abline(h=-bon.logp, v=-bon.logp, col='#cccccc')
    abline(h=0, v=0)
    axis(3, c(-x.thresh, -bon.logp),
         c(x.thresh.lab, expression(p[Bon]==alpha)),
         lwd=0, mgp=c(0, 0, 0), cex.axis=0.6)
    axis(4, c(-y.thresh, -bon.logp),
         c(y.thresh.lab, expression(p[Bon]==alpha)),
         lwd=0, mgp=c(0, 0, 0), cex.axis=0.6)
  }
  if ( add.legend ) {
    legend('topleft', c('Not significant', 'BH sig. in both',
                        paste('BH sig.', x.lab),
                        paste('BH sig.', y.lab)), col=colors, pch=par('pch'), bty='n')
  }
  if ( !missing(filename) ) {
    dev.off()
  }
}
