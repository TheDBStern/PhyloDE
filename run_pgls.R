library(ape)
library(geiger)
library(nlme)
library(phytools)

outfile <- "PhyloDE_results.csv"
write.csv('Gene_Family,Coef,P-value', file=outfile, row.names=FALSE, col.names=FALSE, append=TRUE)

pgls_loop <- function(treefile){
	gene <- strsplit(treefile,'\\.')[[1]][1]
	dat <- read.csv(paste(gene, ".dat.csv", sep=""), header=T)
	phy <- read.tree(treefile)
	phy <- midpoint.root(phy)
	mod <- gls(State ~ Expression, correlation = corMartins(1, phy = phy), data = dat)
	ans <- anova(mod)
	df <- data.frame(gene, mod$coefficients[2], ans$p-value[2])
	write.table(df, file=outfile, row.names=FALSE, col.names=FALSE, sep=",", append=TRUE)
	}

for(i in list.files(pattern='.tre')){
	pgls_loop(i)
	}
	