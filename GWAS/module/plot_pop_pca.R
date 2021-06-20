#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)

prefix<-args[1]
numpc<-as.numeric(args[2])
probcut<-as.numeric(args[3])
if(probcut<0.5 | probcut>1) {
	stop("Incorrect probability cutoff")
}

label<-basename(prefix)

PCA<-read.table(paste0(prefix,".PCs.txt"), header=TRUE, as.is=TRUE)
RefPop<-subset(PCA, Label!=label)
Samp<-subset(PCA, Label==label)

g<-list()
f<-list()
nplots<-(numpc-1)*numpc/2
k<-1

if(file.exists(paste0(prefix,".Ancestry.txt"))) {
	Pred<-read.table(paste0(prefix,".Ancestry.txt"), header=TRUE, as.is=TRUE)
	Pred<-subset(Pred, Label==label)
	npc<-length(grep("PC_", names(Pred)))
	if(npc<2 || numpc>npc) {
		stop("Incorrect number of PCs: ", npc)
	}
	pops<-gsub("^Prob", "", grep("^Prob\\w+", names(Pred), value=TRUE))
	# Infer sample ancestry from Prob
	Samp$Label<-"Unknown"
	for(pop in pops) {
		Samp$Label<-ifelse(Pred[[paste0("Prob",pop)]]>=probcut, pop, Samp$Label)
	}
	for(ii in seq(1, numpc-1)) {
		for(jj in seq(ii+1, numpc)) {	
			g[[k]]<-local({
				i<-ii
				j<-jj
				ggplot(Samp) + 
				geom_point(data=RefPop, aes(x=RefPop[[paste0("PC_",i)]], y=RefPop[[paste0("PC_",j)]], color=RefPop$Label), alpha=0.2, shape=16) +
				geom_point(aes(x=Samp[[paste0("PC_",i)]], y=Samp[[paste0("PC_",j)]], color=Samp$Label), shape=4) +
				scale_color_brewer(type="qual", palette=6) +
				labs(x=paste0("PC_",i), y=paste0("PC_",j), color="Predicted Ancestry")
				})
			f[[k]]<-local({
				i<-ii
				j<-jj
				ggplot(RefPop) +
				geom_point(aes(x=RefPop[[paste0("PC_",i)]], y=RefPop[[paste0("PC_",j)]], color=RefPop$Label), shape=1) +
				scale_color_brewer(type="qual", palette=6) +
				labs(x=paste0("PC_",i), y=paste0("PC_",j), color="Population Ancestry")
				})
			k<-k+1
		}
	}
} else {
	for(ii in seq(1, numpc-1)) {
		for(jj in seq(ii+1, numpc)) {	
			g[[k]]<-local({
				i<-ii
				j<-jj
				ggplot(Samp) + 
				geom_point(data=RefPop, aes(x=RefPop[[paste0("PC_",i)]], y=RefPop[[paste0("PC_",j)]], color=RefPop$Label), shape=2) +
				geom_point(aes(x=Samp[[paste0("PC_",i)]], y=Samp[[paste0("PC_",j)]]), color="grey50", alpha=0.6, shape=4) +
				scale_color_brewer(type="qual", palette=6) +
				labs(x=paste0("PC_",i), y=paste0("PC_",j), color="Reference Ancestry")
				})
			f[[k]]<-local({
				i<-ii
				j<-jj
				ggplot(RefPop) +
				geom_point(aes(x=RefPop[[paste0("PC_",i)]], y=RefPop[[paste0("PC_",j)]], color=RefPop$Label), shape=1) +
				scale_color_brewer(type="qual", palette=6) +
				labs(x=paste0("PC_",i), y=paste0("PC_",j), color="Population Ancestry")
				})
			k<-k+1
		}
	}

}


# Draw pairwise scatterplots

ncol<-ceiling(sqrt(nplots))
nrow<-ceiling(nplots/ncol)

if(ncol>3) {
	width<-15
} else {
	width<-6+(ncol-1)*3
}
height<-nrow/ncol*width
if(height!=width) {
	height <- height+1
}


png(paste0(prefix, '.PCs1-', numpc, "_Prob", probcut, ".png"), width=width, height=height, res=300, unit="in")
ggarrange(plotlist=g, ncol=ncol, nrow=nrow, common.legend = TRUE, legend="bottom")
dev.off()

png(paste0(prefix, '_PopRefs.PCs1-', numpc, ".png"), width=width, height=height, res=300, unit="in")
ggarrange(plotlist=f, ncol=ncol, nrow=nrow, common.legend = TRUE, legend="bottom")
dev.off()

