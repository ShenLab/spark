library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

PREFIX=args[1]

X <- read.table(paste0(PREFIX, "_merged.txt"), header=T, na.string=".")

X$NDP_ChrY <- X$DP_ChrY/X$DP_Auto
X$NDP_ChrX <- X$DP_ChrX/X$DP_Auto
X$NDP_Chr21<- X$DP_Chr21/X$DP_Auto

if(all(is.na(X$NDP_ChrY))) {
  sexchk <- ggplot(X, aes(x=NDP_ChrX, y=0.5, color=Gender)) + geom_jitter(shape=4, size=2)
} else {
  sexchk <- ggplot(X, aes(x=NDP_ChrX, y=NDP_ChrY, color=Gender)) + geom_point(shape=4, size=2)
}

ggsave(paste0(PREFIX, "_chrXY.png"), sexchk, height=5, width=6)

#Err<-subset(X, NDPChrY<0.1 & (NDPChrX<0.8 | NDPChrX>1.2) | NDPChrY>0.1 & NDPChrX>0.8)
#write.csv(Err, "pVCFs/DeepVar/SiteDP_chrXYabnorm.csv")

trisomy <- ggplot(X, aes(NDP_Chr21)) + stat_ecdf(geom="point", shape=4) + ylab("CDF")
ggsave(paste0(PREFIX, "_chr21.png"), trisomy, height=5, width=5)
#Trisomy<-subset(X, NDPChr21>1.2)
#write.csv(Trisomy, "pVCFs/DeepVar/SiteDP_chr21trisomy.csv")


