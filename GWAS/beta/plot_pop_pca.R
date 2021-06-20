#!/usr/bin/env Rscript

library(GGally)

args = commandArgs(trailingOnly=TRUE)

prefix<-args[1]
numpc<-as.numeric(args[2])
testpop<-args[3]


if (numpc < 3) {
	stop("At least 3 PCs should be plotted", call.=FALSE)
}

PCA<-read.table(paste0(prefix,".PCs.txt"), header=TRUE, as.is=TRUE)

Ref<-subset(PCA, Label != testpop)

JJ<-3+as.numeric(numpc)


f<-ggpairs(data=Ref, columns=4:JJ, 
			mapping=ggplot2::aes(color=Label), 
			upper="blank", 
			lower=list(continuous = wrap("points", shape=4, size=0.5)), 
			diag=list(continuous= "barDiag", binwidth=0.5), legend=5) + theme(legend.position="bottom")

width<-as.integer(sqrt(as.numeric(numpc-1))*4)
height<-width+1

ggsave(paste0(prefix, '_PopRefs.PCs1-', numpc, ".png"), f, width = width, height = height, unit = "in")



g<-ggpairs(data=PCA, columns=4:JJ, 
			mapping=ggplot2::aes(color=Label), 
			upper="blank", 
			lower=list(continuous = wrap("points", shape=4, size=0.5)), 
			diag=list(continuous= "barDiag", binwidth=0.5), legend=5) + theme(legend.position="bottom")

width<-as.integer(sqrt(as.numeric(numpc))*3.5)
height<-width+1

ggsave(paste0(prefix, '_Merged.PCs1-', numpc, ".png"), g, width = width, height = height, unit = "in")

