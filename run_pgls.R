# code to run pgls on gene trees (state vs expression) taking into account within-species variance
# for some reason I could not get it to run as a function, so it is a loop

library(ape)
library(geiger)
library(nlme)
library(phytools)
library(stats)

outfile <- "PhyloDE_results.csv"

dfout <- data.frame()
for(treefile in list.files(pattern='.tre')){
	gene <- strsplit(treefile,'\\.')[[1]][1]
	dat <- read.csv(paste(gene, ".dat.csv", sep=""), header=T, row.names=1)
	phy <- read.tree(treefile)
	# adjust branch lengths to avoid computational errors
	phy$edge.length <- phy$edge.length * 100
	phy <- midpoint.root(phy)
	# weighted gls because tree not ultrametric
	w <- diag(vcv.phylo(phy))
	# take into account within-species variance
	vars <- structure(dat[[3]], names=row.names(dat))
	vf <- varComb(varFixed(~ w), varFixed(~ vars))
	mod <- gls(ExpMean ~ State, correlation = corMartins(1, phy = phy, fixed=FALSE), data = dat, weights=vf)
	ans <- anova(mod)
	df <- data.frame(gene, mod$coefficients[[2]], ans$'p-value'[[2]])
	names(df) <- c("Gene","Coefficient","P.Value")
	dfout <- rbind(dfout, df)
	}

# add FDR correction
p <- dfout$P.Value
FDR <- p.adjust(p, "fdr")
dfout <- cbind(dfout, FDR)

write.table(dfout, file=outfile, row.names=FALSE, sep=",")
