#!/usr/bin/env Rscript
library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)

prefix<-args[1]
probcut<-as.numeric(args[2])
if(probcut<0.5 | probcut>1) {
	stop("Incorrect probability cutoff")
}

label<-basename(prefix)

PCA<-read.table(paste0(prefix,".PCs.txt"), header=TRUE, as.is=TRUE)
RefPop<-subset(PCA, Label!=label)


Samp<-read.table(paste0(prefix,".Ancestry.txt"), header=TRUE, as.is=TRUE)

# Determine number of PCs and number populations
npc<-length(grep("PC_", names(Samp)))
if(npc<2) {
	stop("Incorrect number of PCs")
}

pops<-gsub("^Prob", "", grep("^Prob\\w+", names(Samp), value=TRUE))

# Infer sample ancestry from Prob
Samp$Label<-"Unknown"
for(pop in pops) {
	Samp$Label<-ifelse(Samp[[paste0("Prob",pop)]]>=probcut, pop, Samp$Label)
}

nplots<-(npc-1)*npc/2

# Draw pairwise scatterplots
g<-list()
k<-1
for(ii in seq(1, npc-1)) {
	for(jj in seq(ii+1, npc)) {	
		g[[k]]<-local({
			i<-ii
			j<-jj
			ggplot(Samp) + 
				geom_point(aes(x=Samp[[paste0("PC_",i)]], y=Samp[[paste0("PC_",j)]], color=Samp$Label), alpha=0.8, shape=4) +
				geom_point(data=RefPop, aes(x=RefPop[[paste0("PC_",i)]], y=RefPop[[paste0("PC_",j)]], color=RefPop$Label), shape=1, stroke=1, size=4) +
				scale_color_brewer(type="qual", palette=6) +
				labs(x=paste0("PC_",i), y=paste0("PC_",j), color="Predicted Ancestry")
				})
		k<-k+1
	}
}

# Save the discrete population label from prediction
write.table(Samp, file=paste0(prefix,".PopGroup_P", probcut,".txt"), sep="\t", quote=FALSE, row.names=FALSE)

ncol<-ceiling(sqrt(nplots))
nrow<-ceiling(nplots/ncol)

if(ncol>3) {
	width<-15
} else {
	width<-6+(ncol-1)*3
}
height<-ncol/nrow*width
if(height!=width) {
	height <- height+1
}


png(paste0(prefix, '.PCs1-', npc, ".png"), width=width, height=height, res=300, unit="in")
ggarrange(plotlist=g, ncol=ncol, nrow=nrow, common.legend = TRUE, legend="bottom")
dev.off()
