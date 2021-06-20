library(ggplot2)
library(plotrix)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Must provide WRKDIR and PARAMS", call.=FALSE)
}

WRKDIR=args[1]
if (!file_test("-d", WRKDIR)) {
        stop(paste(WRKDIR," is not a directory"), call.=FALSE)
}
PARAMS=args[2]
if (!file_test("-f", PARAMS)) {
        stop(paste(PARAMS, " is not a file"), call.=FALSE)
}

source(PARAMS)

SNP <- read.table(paste0(WRKDIR,"/pass_3.snp_qc_stats.txt"), header=T)

RARE <- subset(SNP, !is.na(MAF) & MAF<=RAREMAF)
COMM <- subset(SNP, !is.na(MAF) & MAF>RAREMAF)


# CDF of SNP call rate separating common and rare SNPs
RARE$CallRate <- 1-RARE$F_MISS
RARE <- RARE[order(RARE$CallRate),]
RARE$Quantile <- (1:dim(RARE)[1])/(dim(RARE)[1]-1)
g <- ggplot(RARE) + geom_point(aes(x=Quantile,y=CallRate)) +
	labs(y="Call rate", title="Ordered coverage for rare SNPs") +
	geom_hline(yintercept=1-MAXRMIS, col="orange", linetype="dashed")
ggsave(file=paste0(WRKDIR,"/pass_3_rare_miss.png"), g, unit="in", width=5, height=5)

COMM$CallRate <- 1-COMM$F_MISS
COMM <- COMM[order(COMM$CallRate),]
COMM$Quantile <- (1:dim(COMM)[1])/(dim(COMM)[1]-1)
g <- ggplot(COMM) + geom_point(aes(x=Quantile,y=CallRate)) +
	labs(y="Call rate", title="Ordered coverage for common SNPs") +
	geom_hline(yintercept=1-MAXCMIS, col="orange", linetype="dashed")
ggsave(file=paste0(WRKDIR,"/pass_3_comm_miss.png"), g, unit="in", width=5, height=5)


# HWE p-value QQ plot
if(CTRLHWE == TRUE) {
	HWEP <- SNP$CTRL_HWE
} else {
	HWEP <- SNP$HWE
}
HWEP <- HWEP[is.finite(HWEP)]
n <- length(HWEP)
png(file=paste0(WRKDIR,"/pass_3_hwe_qq.png"), width=5, height=6, unit="in", res=300)
plot( qchisq((1:n)/(n+1),2), sort(-2*log(HWEP)),
       main="Q-Q plot of -2*ln(HWE P-values)",
       xlab="Expected quantile", ylab="Observed quantile" )
grid()
abline(a=0, b=1, col="red")
abline(h=-2*log(MINHWE), col="orange", lty=3)
dev.off()


gapped_barplot <- function(VAR, Thres, File, Atten=0.1) {
	height <- table(SNP[[VAR]])
	if(length(height) <= 2) {
		png(file=File, width=6, height=6, unit="in", res=300)
		plot(0,type='n',axes=FALSE,ann=FALSE)
		dev.off()
	}
	htsort <- sort(height, decreasing=TRUE)
	scl  <- (htsort[1]/htsort[2])^Atten
	from <- htsort[2]*scl
	to   <- htsort[1]/scl
	png(file=File, width=6, height=6, unit="in", res=300)
	gap.barplot(height, gap=c(from, to), ylim=c(0,2*from), ytics=c(as.integer(seq(1,2)*from/2), as.integer(to+seq(1,2)*from/2)), xlab=VAR, ylab="Counts")
	title(paste("Histogram of",VAR))
	axis.break(2, from, breakcol="snow", style="gap")
	axis.break(2, from*(1+0.02), breakcol="black", style="slash")
	axis.break(4, from*(1+0.02), breakcol="black", style="slash")
	abline(v=Thres, col="orange")
	dev.off()
}

# Histogram of HH count
if(file.exists(paste0(WRKDIR,"/pass_3_hist_hhcount.png")))
	file.remove(paste0(WRKDIR,"/pass_3_hist_hhcount.png"))
if ("HH_COUNT" %in% names(SNP)) {
	gapped_barplot("HH_COUNT", MAXHHCNT, paste0(WRKDIR,"/pass_3_hist_hhcount.png"), 0.1)
} else {
	png(paste0(WRKDIR,"/pass_3_hist_hhcount.png"), width=5, height=5, unit="in", res=300)
	plot(0,type='n',axes=FALSE,ann=FALSE)
	dev.off()
}

# Histogram of Mendel error
if(file.exists(paste0(WRKDIR,"/pass_3_hist_mendel.png")))
	file.remove(paste0(WRKDIR,"/pass_3_hist_mendel.png"))
if ("MENDEL" %in% names(SNP)) {
	gapped_barplot("MENDEL", MAXMENDEL, paste0(WRKDIR,"/pass_3_hist_mendel.png"), 0.02)
}

# Histogram of duplication error
if(file.exists(paste0(WRKDIR,"/pass_3_hist_duperr.png")))
	file.remove(paste0(WRKDIR,"/pass_3_hist_duperr.png"))
if ("DUPERR" %in% names(SNP)) {
	gapped_barplot("DUPERR", MAXDUPERR, paste0(WRKDIR,"/pass_3_hist_duperr.png"), 0.01)
} else {
	png(paste0(WRKDIR,"/pass_3_hist_duperr.png"), width=5, height=5, unit="in", res=300)
	plot(0,type='n',axes=FALSE,ann=FALSE)
	dev.off()
}

# Now create SNP exclusion list
MASK <- with(SNP, is.na(MAF) | MAF<MINMAF | MAF<=RAREMAF & F_MISS>MAXRMIS |
			 	  MAF>RAREMAF & F_MISS>MAXCMIS | !is.na(HWE) & HWE<MINHWE)
if("HH_COUNT" %in% names(SNP)) {
	EXCL <- with(SNP, !is.na(HH_COUNT) & HH_COUNT>MAXHHCNT)
	MASK <- MASK | EXCL
}
if("MENDEL" %in% names(SNP)) {
	EXCL <- with(SNP, !is.na(MENDEL) & MENDEL>MAXMENDEL)
	MASK <- MASK | EXCL
}
if("DUPERR" %in% names(SNP)) {
	HQ <- with(SNP, !is.na(DUPERR) & DUPERR>MAXDUPERR)
	MASK <- MASK | EXCL
}

BADSNP <- SNP[MASK,]
write.table(BADSNP, paste0(WRKDIR, "/pass_3.badsnps.txt"), sep="\t", quote=FALSE, row.names=FALSE)

EXCLUDE <- BADSNP$SNP
write.table(EXCLUDE, paste0(WRKDIR, "/pass_3.exclude"), row.names=FALSE, quote=FALSE, col.names=FALSE)

