# Pick columns (phenotypes) at random and make a bunch of QQ plots

if ( length(commandArgs(TRUE)) < 1 & !exists('PREFIX')) {
  stop('Please specify a data file on the command line or set the PREFIX variable')
} else if ( !exists('PREFIX') ) {
  PREFIX <<- commandArgs(TRUE)[1]
}

require('bigmemory')

data.file <- sprintf('%s_pvals.desc', PREFIX)

if ( !exists('pvals') ) {
	cat('Loading', data.file, '...\n')
	pvals <- attach.big.matrix(data.file)
}

if ( is.null(colnames(pvals)) ) {
	columns <- sample(1:ncol(pvals), min(80, ncol(pvals)))
} else {
	columns <- sample(colnames(pvals), min(80, ncol(pvals)))
}

source('qqplots.r')

#make.qq.plots(columns, pvals, paste(PREFIX, '_random', sep=''))


library("bigmemory")

pvalues_matrix<-attach.big.matrix(data.file)

#Q-Q_plot function for the analysis

ggd.qqplot = function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,max(o)))
    lines(e,e,col="red")
}


png("P_values_diagnostic_ALL_TE_INSERTION_ALL_EXPRESSION_plot.png")

par(cex.axis=1.25, lwd=3)
ggd.qqplot(as.numeric(as.vector(as.matrix(pvalues_matrix))))
axis(side = 2, font = 2)
axis(side = 1, font = 2)
dev.off()


save.image("QQ_plot.Rdat")