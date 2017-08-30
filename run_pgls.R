library(ape)
library(geiger)
library(nlme)
library(phytools)
library(stats)

outfile <- "PhyloDE_results.csv"
writeLines('Gene,Coefficient,P-Value', outfile)

pgls_loop <- function(treefile){
	gene <- strsplit(treefile,'\\.')[[1]][1]
	dat <- read.csv(paste(gene, ".dat.csv", sep=""), header=T, row.names=1)
	phy <- read.tree(treefile)
	# adjust branch lengths to avoid computational errors
	phy$edge.length <- phy$edge.length * 100
	phy <- midpoint.root(phy)
	# weighted gls because tree not ultrametric
	w<-diag(vcv.phylo(phy))
	# take into account within-species variance
	vars <- structure(dat[[3]], names=row.names(dat))
	vf <- varComb(varFixed(~ w), varFixed(~ vars))
	mod <- gls(ExpMean ~ State, correlation = corMartins(1, phy = phy, fixed=FALSE), data = dat, weights=vf)
	ans <- anova(mod)
	df <- data.frame(gene, mod$coefficients[2], ans$'p-value'[2])
	write.table(df, file=outfile, row.names=FALSE, sep=",", append=TRUE)
	}

for(i in list.files(pattern='.tre')){
	pgls_loop(i)
	}


