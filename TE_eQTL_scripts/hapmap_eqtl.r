if ( length(commandArgs(TRUE)) < 1 & !exists('pop')) {
  stop('Please specify a population identifier in the command line...')
} else if ( !exists('pop') ) {
  pop <- commandArgs(TRUE)[1]
}


data.file <- 'hapmap_eqtl.RData'
load(data.file)

all.genotypes <- genotypes

if ( ! pop %in% pops ) {
  cat('Found pops:', pops, '\n')
  stop(pop, ' was not found in the data file.\n')
}

phenotypes <- t(expression[, samples$id[samples$pop==pop]])
n <- pop.counts[pop]
covars <- samples[samples$id[samples$pop==pop], 'sex']
genotypes <- t(all.genotypes[use.snps[, pop], samples$id[samples$pop==pop]])



PREFIX <- sprintf('%s_hapmap', pop)
MC.CORES <- 24
source('mc_wald_xqtl.r')