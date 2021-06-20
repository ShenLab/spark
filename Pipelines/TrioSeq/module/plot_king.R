library(ggplot2)
#library(ggpubr)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

PREFIX=args[1]

X<-read.table(paste0(PREFIX, ".kin.txt"), as.is=T, header=T)
Y<-read.table(paste0(PREFIX, ".kin0.txt"), as.is=T, header=T)

# Collect cyptic relpairs
Cyp<-subset(Y, InfType=="Dup/MZ" | InfType=="FS" | InfType=="PO")
Cyp$PedType<-Cyp$InfType

#table(X$Relpair)

g1<-ggplot(X, aes(x=IBS0, y=Kinship, color=PedType)) + geom_point(shape="x", size=2) +
              geom_point(data=Cyp, aes(x=IBS0, y=Kinship, color=PedType)) +
              theme(legend.position="none")
g2<-ggplot(X, aes(x=IBD2Seg, y=Kinship, color=PedType)) + geom_point(shape="x", size=2) +
              geom_point(data=Cyp, aes(x=IBD2Seg, y=Kinship, color=PedType)) 

png(paste0(PREFIX, ".png"), height=5, width=11, unit="in", res=300)
#ggarrange(g1, g2, ncol=2, common.legend=TRUE, legend="top")
grid.arrange(g1, g2, ncol=2, widths=c(4, 5))
dev.off()

#g3<-ggplot(X, aes(x=IBS0, y=IBD1Seg, color=PedType)) + geom_point(shape="x", size=2) +
#              geom_point(data=Cyp, aes(x=IBS0, y=IBD1Seg, color=PedType)) 
#g4<-ggplot(X, aes(x=IBD2Seg, y=IBD1Seg, color=PedType)) + geom_point(shape="x", size=2) +
#              geom_point(data=Cyp, aes(x=IBD2Seg, y=IBD1Seg, color=PedType)) 

#png(paste0(PREFIX, ".png"), height=12, width=12, unit="in", res=300)
#ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, common.legend=TRUE, legend="right")
#dev.off()


