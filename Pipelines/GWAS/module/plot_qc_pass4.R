library(ggplot2)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3 && length(args) !=4) {
  stop("Must provide PREFIX and PARAMS", call.=FALSE)
}

PREFIX=args[1]
PEDFILE=args[2]
PARAMS=args[3]
if (!file_test("-f", PARAMS)) {
	stop(paste(PARAMS, " is not a file"), call.=FALSE)
}
if (length(args) == 3) {
	OUTFILE=PREFIX
} else {
	OUTFILE=args[4]
}

source(PARAMS)

PED<-subset(read.table(PEDFILE, header=F, as.is=T),select=c("V1","V2","V5"));
names(PED)<-c("FID","IID","PEDSEX")
PED$PEDSEX<-with(PED, ifelse(PEDSEX==1, "Male", ifelse(PEDSEX==2, "Female", "Unknown")))


IMISS<-subset(read.table(paste0(PREFIX, ".imiss"), header=T, as.is=T),
	 			select=c("FID","IID","F_MISS"))
names(IMISS)[3]<-"GMISS"
HET<-subset(read.table(paste0(PREFIX, ".het"), header=T, as.is=T),
				select=c("FID","IID","F"))
names(HET)[3]<-"F_AUTO"

if (file.exists(paste0(PREFIX, "_chrXY.sexcheck"))) {
	SEX<-subset(read.table(paste0(PREFIX, "_chrXY.sexcheck"), header=T, as.is=T),
	 			select=c("FID","IID","F","YCOUNT"))
} else {
	SEX<-subset(read.table(paste0(PREFIX, "_chrX.sexcheck"), header=T, as.is=T),
	 			select=c("FID","IID","F"))
}
names(SEX)[3]<-"F_X"
SEX$SNPSEX<-with(SEX, ifelse(F_X<=FFMAX, "Female", ifelse(F_X>=MFMIN, "Male", "Unknown")))

SUM<-merge(merge(merge(IMISS, HET), PED), SEX)
write.csv(SUM, file=paste0(OUTFILE, "_sample.csv"), quote=FALSE, row.names=FALSE)

# Bad sample: high missing rate or unusual F 
BADSAMP<-subset(SUM, GMISS>MISSMAX | F_AUTO<WGFMIN | F_AUTO>WGFMAX)

if(file.exists(paste0(OUTFILE, ".badsamp.csv")))
	file.remove(paste0(OUTFILE, ".badsamp.csv"))
if(file.exists(paste0(OUTFILE, ".badsamp.remove")))
	file.remove(paste0(OUTFILE, ".badsamp.remove"))
if (nrow(BADSAMP) > 0) {
	write.csv(BADSAMP, file=paste0(OUTFILE, ".badsamp.csv"), quote=FALSE, row.names=FALSE)
	write.table(BADSAMP[,c("FID","IID")], file=paste0(OUTFILE, ".badsamp.remove"), quote=FALSE, 
		row.names=FALSE, col.names=FALSE)
}

#
# CDF of call rate and inbreeding coeff
#
g<-ggplot(SUM) + geom_point(aes(x=F_AUTO, y=GMISS)) +
	geom_hline(yintercept=MISSMAX, color="orange", linetype="dashed") +
	geom_vline(xintercept=WGFMIN, color="darkgreen", linetype="dashed") +
	geom_vline(xintercept=WGFMAX, color="darkgreen", linetype="dashed") +
	labs(x="Inbreeding coeff", y="Missing genotype rate")
if (nrow(BADSAMP)>0) {
	g<-g+geom_point(data=BADSAMP, aes(x=F_AUTO, y=GMISS))+
		geom_text_repel(data=BADSAMP, aes(x=F_AUTO, y=GMISS), label=BADSAMP$IID, size=3)
}
ggsave(paste0(OUTFILE, "_HetvsMiss.png"), g, unit="in", width=6, height=5)

# Update Sex:
# Update sex based on chrX het rate, unknown SNPSEX will be 0
SEXUPDATE<-subset(SUM, PEDSEX!=SNPSEX & ! paste(FID, IID) %in% rownames(BADSAMP), 
		select=c("FID","IID","SNPSEX"))
if(file.exists(paste0(OUTFILE, ".sexupdate")))
	file.remove(paste0(OUTFILE, ".sexupdate"))
if (nrow(SEXUPDATE) > 0) {
	write.table(SEXUPDATE, paste0(OUTFILE, ".sexupdate"), quote=FALSE,
		col.names=FALSE, row.names=FALSE)
}

SEXERR<-subset(SUM, PEDSEX!=SNPSEX)
if(file.exists(paste0(OUTFILE, ".sexerr.csv")))
	file.remove(paste0(OUTFILE, ".sexerr.csv"))
if (nrow(SEXERR) > 0) {
	write.csv(SEXERR, file=paste0(OUTFILE, ".sexerr.csv"), quote=FALSE, row.names=FALSE)
}

# Remove sample with unknown sex?
if(file.exists(paste0(OUTFILE, ".sexna.remove")))
	file.remove(paste0(OUTFILE, ".sexna.remove"))
if (exists("SEXIMP") && SEXIMP == TRUE) {
	SEXERR<-subset(SUM, SNPSEX=="Unknown")
	write.table(SEXERR[,c("FID","IID")], file=paste0(OUTFILE, ".sexerr.remove"), quote=FALSE, 
		row.names=FALSE, col.names=FALSE)
}


#
# chrX F vs whole genome missing rate
#
g<-ggplot(SUM) + geom_point(aes(x=GMISS, y=F_X, color=PEDSEX, shape=SNPSEX)) +
	geom_hline(yintercept=FFMAX, color="darkgreen", linetype="dashed") +
	geom_hline(yintercept=MFMIN, color="darkgreen", linetype="dashed") +
	labs(x="Missing genotype rate", y="chrX inbreeding coeff", color="PedSex", shape="SNPSex")
if (nrow(SEXERR) > 0) {
	g<-g+geom_point(data=SEXERR, aes(x=GMISS, y=F_X, color=PEDSEX, shape=SNPSEX))+
		geom_text_repel(data=SEXERR, aes(x=GMISS, y=F_X), label=SEXERR$IID, size=3)
}
ggsave(paste0(OUTFILE, "_chrX_FXvsMiss.png"), g, unit="in", width=6, height=5)

#
# chrX F vs autosome F.  
#
g<-ggplot(SUM) + geom_point(aes(x=F_AUTO, y=F_X, color=PEDSEX, shape=SNPSEX)) +
	geom_hline(yintercept=FFMAX, color="darkgreen", linetype="dashed") +
	geom_hline(yintercept=MFMIN, color="darkgreen", linetype="dashed") +
	labs(x="Autosome inbreeding coeff", y="chrX inbreeding coeff", color="PedSex", shape="SNPSex")
if (nrow(SEXERR) > 0) {
	g<-g+geom_point(data=SEXERR, aes(x=F_AUTO, y=F_X, color=PEDSEX, shape=SNPSEX))+
		geom_text_repel(data=SEXERR, aes(x=F_AUTO, y=F_X), label=SEXERR$IID, size=3)
}
ggsave(paste0(OUTFILE, "_chrX_F_XvsAUTO.png"), g, unit="in", width=6, height=5)


# If further chrY markers are available
#
# chrX F vs Y geno counts
#
if(file.exists(paste0(OUTFILE, "_chrXY_FXvsYcnt.png")))
	file.remove(paste0(OUTFILE, "_chrXY_FXvsYcnt.png"))
if("YCOUNT" %in% names(SUM)) {
	g<-ggplot(SUM) + geom_point(aes(x=YCOUNT, y=F_X, color=PEDSEX, shape=PEDSEX)) +
		geom_hline(yintercept=FFMAX, color="darkgreen", linetype="dashed") +
		geom_hline(yintercept=MFMIN, color="darkgreen", linetype="dashed") +
		labs(x="chrY genotype counts", y="chrX inbreeding coeff", color="PedSex", shape="SNPSex")
	if (nrow(SEXERR) > 0) {
		g<-g+geom_point(data=SEXERR, aes(x=YCOUNT, y=F_X, color=PEDSEX, shape=SNPSEX))+
			geom_text_repel(data=SEXERR, aes(x=YCOUNT, y=F_X), label=SEXERR$IID, size=3)
	}
	ggsave(paste0(OUTFILE, "_chrXY_FXvsYcnt.png"), g, unit="in", width=6, height=5)
}

# Write final list of bad samples

