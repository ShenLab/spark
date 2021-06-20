library(ggplot2)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2 && length(args)!=3) {
  stop("Must provide PREFIX and PARAMS", call.=FALSE)
}

PREFIX=args[1]
PARAMS=args[2]
if (!file_test("-f", PARAMS)) {
	stop(paste(PARAMS, " is not a file"), call.=FALSE)
}
if (length(args) == 2) {
	OUTFILE=PREFIX
} else {
	OUTFILE=args[3]
}

source(PARAMS)

# Plot CDF of missing call rate
IMISS<-read.table(paste0(PREFIX, ".imiss"), header=T, as.is=T)
IMISS$CALLRATE<-1-IMISS$F_MISS
IMISS<-IMISS[order(IMISS$CALLRATE),]
IMISS$QUANTILE<-(1:dim(IMISS)[1])/(dim(IMISS)[1]-1)

BADSAMP<-subset(IMISS, F_MISS>MISSMAX)

#
# CDF of call rate
#
g<-ggplot(IMISS) + geom_point(aes(x=QUANTILE, y=CALLRATE)) +
	geom_hline(yintercept=1-MISSMAX, col="orange", linetype="dashed") +
	labs(x="Quantile", y="Call rate", title="Ordered individual call rate")
if(file.exists(paste0(OUTFILE, "_miss.badsamps.csv")))
	file.remove(paste0(OUTFILE, "_miss.badsamps.csv"))
if(nrow(BADSAMP)>0) {
	# Avoid displaying too many bad samples
	if (nrow(BADSAMP)<20) {
		g<-g+geom_text_repel(data=BADSAMP, aes(x=QUANTILE, y=CALLRATE), label=BADSAMP$IID, size=3)
	}
	write.csv(BADSAMP, file=paste0(OUTFILE, "_miss.badsamps.csv"), quote=FALSE, row.names=FALSE)
}
ggsave(file=paste0(OUTFILE, "_miss.png"), g, unit="in", width=6, height=5)

# Scatterplot of inbreeding coeff vs missing rate using pre-QC SNPs
# Can be used to spot sample contamination
HET<-read.table(paste0(PREFIX, ".het"), header=T, as.is=T)
HETMIS<-merge(HET, IMISS)

# Check for samples with large het or low F
CONTSAMP<-subset(HETMIS, F < WGFMIN)

#
# Inbreeding vs missing rate
#
g<-ggplot(HETMIS) + geom_point(aes(x=F, y=F_MISS)) +
	geom_vline(xintercept=WGFMIN, col="orange", linetype="dashed") +
	geom_hline(yintercept=MISSMAX, col="darkgreen", linetype="dashed") +
	labs(x="Inbreeding coeff estimate", y="Missing genotype rate")
if(file.exists(paste0(OUTFILE, "_het.badsamps.csv")))
	file.remove(paste0(OUTFILE, "_het.badsamps.csv"))
if (nrow(CONTSAMP) > 0) {
	g<-g+geom_text_repel(data=CONTSAMP, aes(x=F, y=F_MISS), label=CONTSAMP$IID, size=3)
	write.csv(CONTSAMP, file=paste0(OUTFILE, "_het.badsamps.csv"), quote=FALSE, row.names=FALSE)
} 
if (nrow(BADSAMP) > 0) {
	BADSAMP2<-subset(HETMIS, F>=WGFMIN & F_MISS>MISSMAX)
	if (nrow(BADSAMP2) > 0) {
		g<-g+geom_text_repel(data=BADSAMP2, aes(x=F, y=F_MISS), label=BADSAMP2$IID, size=3)
	}
}
ggsave(file=paste0(OUTFILE,"_het.png"), g, unit="in", width=6, height=5)


# Gender vs chrX het, also show missing rate
SEX<-read.table(paste0(PREFIX, "_chrX.sexcheck"), header=T)
SEX$SNPSex<-ifelse(SEX$F <= FFMAX, "Female", ifelse(SEX$F >= MFMIN, "Male", "Unknown"))
SEX$SNPSEX<-ifelse(SEX$F <= FFMAX, 2, ifelse(SEX$F >= MFMIN, 1, 0))
SEX$PedSex<-ifelse(SEX$PEDSEX == 1, "Male", ifelse(SEX$PEDSEX == 2, "Female", "Unknown"))


SEXMIS<-merge(SEX, IMISS)

# Check for problematic samples
SEXERR<-subset(SEXMIS, PedSex!="Unknown" & SNPSex!=PedSex)

#
# chrX F and missing rate, and highlight sample with problematic sex
#
g<-ggplot(SEXMIS) + geom_point(aes(x=F_MISS, y=F, color=PedSex, shape=SNPSex)) +
	geom_vline(xintercept=MISSMAX, color="orange", linetype="dashed") +
	geom_hline(yintercept=FFMAX, color="darkgreen", linetype="dashed") +
	geom_hline(yintercept=MFMIN, color="darkgreen", linetype="dashed") +
	labs(x="Missing genotype rate", y="chrX inbreeding coeff")
if(file.exists(paste0(OUTFILE, "_chrX.sexerr.csv")))
	file.remove(paste0(OUTFILE, "_chrX.sexerr.csv"))
if (nrow(SEXERR) > 0) {
	g<-g+geom_point(data=SEXERR, aes(x=F_MISS, y=F, color=PedSex, shape=SNPSex))+
		geom_text_repel(data=SEXERR, aes(x=F_MISS, y=F), label=SEXERR$IID, size=3)
	write.csv(SEXERR, file=paste0(OUTFILE, "_chrX.sexerr.csv"), quote=FALSE, row.names=FALSE)
}
ggsave(paste0(OUTFILE, "_chrX_FvsMiss.png"), g, unit="in", width=6, height=5)


# If further chrY markers are available

#
# chrY geno counts vs chrX F
#
if(file.exists(paste0(OUTFILE, "_chrXY_FXvsYcnt.png")))
	file.remove(paste0(OUTFILE, "_chrXY_FXvsYcnt.png"))
if(file_test("-f", paste0(PREFIX, "_chrXY.sexcheck"))) {
	cat("Processing XY", fill=TRUE)
	XY<-read.table(paste0(PREFIX, "_chrXY.sexcheck"), header=T)
	XY$PedSex<-ifelse(XY$PEDSEX == 1, "Male", ifelse(XY$PEDSEX == 2, "Female", "Unknown") )
	XY$SNPSex<-ifelse(XY$F <= FFMAX, "Female", ifelse(XY$F >= MFMIN, "Male", "Unknown"))

	SEXMIS<-merge(XY, IMISS)
	SEXERR<-subset(SEXMIS, PedSex!="Unknown" & SNPSex!=PedSex)

	g<-ggplot(XY) + geom_point(aes(x=YCOUNT, y=F, color=PedSex, shape=SNPSex)) +
		geom_hline(yintercept=FFMAX, color="darkgreen", linetype="dashed") +
		geom_hline(yintercept=MFMIN, color="darkgreen", linetype="dashed") +
		labs(x="chrY genotype counts", y="chrX inbreeding coeff")
	if (nrow(SEXERR) > 0) {
		g<-g+geom_point(data=SEXERR, aes(x=YCOUNT, y=F, color=PedSex, shape=SNPSex))+
			geom_text_repel(data=SEXERR, aes(x=YCOUNT, y=F), label=SEXERR$IID, size=3)
	}
	ggsave(paste0(OUTFILE, "_chrXY_FXvsYcnt.png"), g, unit="in", width=6, height=5)


	g<-ggplot(SEXMIS) + geom_point(aes(x=YCOUNT, y=F_MISS, color=PedSex, shape=SNPSex)) +
		geom_hline(yintercept=MISSMAX, color="darkgreen", linetype="dashed") +
		labs(x="chrY genotype counts", y="Missing genotype rate")
	if (nrow(SEXERR) > 0) {
		g<-g+geom_point(data=SEXERR, aes(x=YCOUNT, y=F_MISS, color=PedSex, shape=SNPSex))+
			geom_text_repel(data=SEXERR, aes(x=YCOUNT, y=F_MISS), label=SEXERR$IID, size=3)
	}
	ggsave(paste0(OUTFILE, "_chrXY_MissvsYcnt.png"), g, unit="in", width=6, height=5)
}

#SEXNA<-subset(SEX, SNPSex == "Unknown")


# Create a full list of bad samples for removal
# They include: samples with low call rate, high heterozygosity
# Note: sex unknown case will not be removed
REMOVE<-unique(rbind(BADSAMP[c("FID","IID")], CONTSAMP[c("FID","IID")]))
if(file.exists(paste0(OUTFILE, ".remove")))
	file.remove(paste0(OUTFILE, ".remove"))
if(nrow(REMOVE) > 0) {
	write.table(REMOVE, paste0(OUTFILE, ".remove"), quote=FALSE, 
		col.names=FALSE, row.names=FALSE)
}
rownames(REMOVE)<-paste(REMOVE$FID, REMOVE$IID)

# Update sex based on chrX het rate, unknown SNPSEX will be 0
SEXUPDATE<-subset(SEX, PEDSEX!=SNPSEX & ! paste(FID, IID) %in% rownames(REMOVE), 
		select=c("FID","IID","SNPSEX"))

# If we update sex for parents, we need to fix fam files later
if(file.exists(paste0(OUTFILE, ".sexupdate")))
	file.remove(paste0(OUTFILE, ".sexupdate"))
if (nrow(SEXUPDATE) > 0) {
	write.table(SEXUPDATE, paste0(OUTFILE, ".sexupdate"), quote=FALSE,
		col.names=FALSE, row.names=FALSE)
}

