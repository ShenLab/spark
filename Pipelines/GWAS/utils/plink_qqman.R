library(qqman)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2 && length(args) !=3) {
  stop("Must provide ASSOC BIMFILE and (optionally) OUTFILE", call.=FALSE)
}

ASSOC=args[1]
if (!file_test("-f", ASSOC)) {
	stop(paste(ASSOC, " is not a file"), call.=FALSE)
}

BIMFILE=args[2]
if (!file_test("-f", BIMFILE)) {
	stop(paste(BIMFILE, " is not a file"), call.=FALSE)
}

if (length(args) == 2) {
	OUTFILE=ASSOC
} else {
	OUTFILE=args[3]
}

STAT<-read.table(ASSOC, header=TRUE)
if ("TEST" %in% names(STAT)) {
	STAT<-subset(STAT, TEST=="ADD")
}

BIM<-subset(read.table(BIMFILE, header=FALSE), select=c("V1","V2","V4"))
names(BIM)<-c("CHR","SNP","BP")

X<-merge(STAT, BIM)
cat("Number of markers with P and CHR/BP", dim(X)[1], fill=TRUE)

png(paste0(OUTFILE, ".qq.png"), width=6, height=6, unit="in", res=300)
qq(X$P, main="Q-Q plot of GWAS p-values")
dev.off()

png(paste0(OUTFILE, ".mht.png"), width=12, height=5, unit="in", res=300)
manhattan(X, main="Genome-wide association signals", chrlabs=as.character(sort(unique(X$CHR))))
dev.off()
