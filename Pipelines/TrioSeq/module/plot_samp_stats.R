#

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Must provide OUTDIR and PARAMS", call.=FALSE)
}

OUTDIR=args[1]
if (!file_test("-d", OUTDIR)) {
	stop(paste(OUTDIR," is not a directory"), call.=FALSE)
}
WRKDIR=args[2]
if (!file_test("-d", OUTDIR)) {
	stop(paste(OUTDIR," is not a directory"), call.=FALSE)
}
PARAMS=args[3]
if (!file_test("-f", PARAMS)) {
	stop(paste(PARAMS, " is not a file"), call.=FALSE)
}
source(PARAMS)

SAMP<-read.csv(paste0(OUTDIR,"/sample.csv"), header=T, stringsAsFactors=F)
SAMP$NDP_chrX<-with(SAMP, DP_chrX/DP_auto)
SAMP$NDP_chrY<-with(SAMP, DP_chrY/DP_auto)
SAMP$Pheno<-as.factor(SAMP$Pheno)
SAMP$LOWDP<-with(SAMP, DP_het<MINDP)

# We will determine sex based on chrY normalized depth
if (MINMYDP <= 0 || all(is.na(SAMP$NDP_chrY))) {
	SAMP$SEXERR<-with(SAMP, Sex_ped != Sex_pred)
	SEXERR<-subset(SAMP, Sex_ped != "unknown" & Sex_ped != Sex_pred)
} else {
	SAMP$Sex_pred_Ydp<-ifelse(SAMP$NDP_chrY >= MINMYDP,  "male", "female")
	SAMP$SEXERR<-with(SAMP, Sex_ped != Sex_pred_Ydp)
	SEXERR<-subset(SAMP, Sex_ped != "unknown" & Sex_ped != Sex_pred_Ydp)
}
ERRSAMP<-subset(SAMP, LOWDP | SEXERR)
LOWDP<-subset(SAMP, DP_het < MINDP)

if(file.exists(paste0(OUTDIR,"/sample_err.csv")))
	file.remove(paste0(OUTDIR,"/sample_err.csv"))
if (nrow(ERRSAMP) > 0) {
	write.csv(ERRSAMP, paste0(OUTDIR,"/sample_err.csv"), quote=FALSE, row.names=FALSE)
}


library(ggplot2)
library(grid)
library(ggtern)
library(GGally)
library(ggrepel)


ggsave.png <- function (filename, g, dpi) {
	gt = ggplot_gtable(ggplot_build(g))
	gt$layout$clip[gt$layout$name=="panel"] = "off"
	png(file=filename, height=8, width=10, unit="in", res=300)
	grid.draw(gt)
	dev.off()
	pdf(NULL)
}


# Sample level plots
{
	pdf(NULL)
	g<-ggplot(SAMP) + geom_point(aes(x=GMiss_HQ, y=DP_auto, color=Pheno)) +
	labs(x="Genotype Missingness (HQ SNPs)", y="Mean Autosome Depth", color="Phenotype")
	if(nrow(LOWDP)) {
		#g<-g+geom_text(data=LOWDP, aes(x=GMiss_HQ, y=DP_auto), label=LOWDP$IID, hjust=0, size=2, check_overlap=TRUE)
		g<-g+geom_text_repel(data=LOWDP, aes(x=GMiss_HQ, y=DP_auto), label=LOWDP$IID, size=2)
	}
	ggsave.png(paste0(OUTDIR,"/sample_DPvsGMissHQ.png"), g, dpi=300)


	g<-ggplot(SAMP) + geom_point(aes(x=GMiss_HQ, y=DP_het, color=Pheno)) +
	labs(x="Genotype Missingness (HQ SNPs)", y="Mean Heterozygotes Depth", color="Phenotype")
	if(nrow(LOWDP)) {
		#g<-g+geom_text(data=LOWDP, aes(x=GMiss_HQ, y=DP_het), label=LOWDP$IID, hjust=0, size=2, check_overlap=TRUE)
		g<-g+geom_text_repel(data=LOWDP, aes(x=GMiss_HQ, y=DP_het), label=LOWDP$IID, size=2)
	}
	ggsave.png(paste0(OUTDIR,"/sample_HetDpvsGMissHQ.png"), g, dpi=300)


	g<-ggplot(SAMP) + geom_point(aes(x=NDP_chrX, y=NDP_chrY, color=Sex_ped)) +
	labs(x="Normalized chrX depth", y="Normalized chrY depth", color="Pedigree Sex")
	if(nrow(SEXERR)) {
		#g<-g+geom_text(data=SEXERR, aes(x=NDP_chrX, y=NDP_chrY), label=SEXERR$IID, hjust=0, size=2, check_overlap=TRUE)
		g<-g+geom_text_repel(data=SEXERR, aes(x=NDP_chrX, y=NDP_chrY), label=SEXERR$IID, size=2)
	}
	ggsave.png(paste0(OUTDIR,"/sample_chrXYNormDp_byPEDSex.png"), g, dpi=300)

	g<-ggplot(SAMP) + geom_point(aes(x=NDP_chrX, y=NDP_chrY, color=Sex_pred)) +
		labs(x="Normalized chrX depth", y="Normalized chrY depth", color="Predicted Sex")
	if(nrow(SEXERR)) {
		#g<-g+geom_text(data=SEXERR, aes(x=NDP_chrX, y=NDP_chrY), label=SEXERR$IID, hjust=0, size=2, check_overlap=TRUE)
		g<-g+geom_text_repel(data=SEXERR, aes(x=NDP_chrX, y=NDP_chrY), label=SEXERR$IID, size=2)
	}
	ggsave.png(paste0(OUTDIR,"/sample_chrXYNormDp_byPredSex.png"), g, dpi=300)


	g<-ggplot(SAMP) + geom_point(aes(x=F_chrX, y=HetRt_chrX, color=Sex_pred)) +
		labs(x="chrX Inbreeding Coeff", y="chrX Het/Hom Ratio", color="Predicted Sex")
	if(nrow(SEXERR)) {
		#g<-g+geom_text(data=SEXERR, aes(x=F_chrX, y=HetRt_chrX), label=SEXERR$IID, hjust=0, size=2, check_overlap=TRUE)
		g<-g+geom_text_repel(data=SEXERR, aes(x=F_chrX, y=HetRt_chrX), label=SEXERR$IID, size=2)
	}
	ggsave.png(paste0(OUTDIR,"/sample_chrXFHet_byPredSex.png"), g, dpi=300)

	g<-ggplot(SAMP) + geom_point(aes(x=F_chrX, y=HetRt_chrX, color=Sex_ped)) +
		labs(x="chrX Inbreeding Coeff", y="chrX Het/Hom Ratio", color="Pedigree Sex")
	if(nrow(SEXERR)) {
		#g<-g+geom_text(data=SEXERR, aes(x=F_chrX, y=HetRt_chrX), label=SEXERR$IID, hjust=0, size=2, check_overlap=TRUE)
		g<-g+geom_text_repel(data=SEXERR, aes(x=F_chrX, y=HetRt_chrX), label=SEXERR$IID, size=2)
	}
	ggsave.png(paste0(OUTDIR,"/sample_chrXFHet_byPEDSex.png"), g, dpi=300)

	g<-ggplot(SAMP) + geom_point(aes(x=F_auto, y=HetRt_auto, color=Popu_pred)) +
		labs(x="Autosome Inbreeding Coeff", y="Autosome Het/Hom Ratio", color="Predicted Population")
	ggsave.png(paste0(OUTDIR,"/sample_autoFHet_byPredPop.png"), g, dpi=300)

	#PCA
	# first draw background
	BG<-read.csv(paste0(WRKDIR,"/peddy.background_pca.csv"), header=T)
	g<-ggpairs( data=BG, columns=1:4, 
		  	mapping = ggplot2::aes(color=ancestry), 
		  	upper="blank", 
		  	lower = list(continuous = wrap("points", alpha = 0.1)),
		  	diag=list(continuous= "barDiag", binwidth=0.5, alpha=0.2),
		  	legend = 5)

	grid<-subset(expand.grid(x=2:4,y=1:3), x>y)

	# then superimpose current sample
	SAMP$ancestry<-SAMP$Popu_pred
	for(ii in 1:nrow(grid)) {
		jj<-grid[ii,1]
		kk<-grid[ii,2]
		cat(jj,kk,fill=TRUE)
		g[jj,kk]<-g[jj,kk]+
			geom_point(data=SAMP, 
				aes_string(x=paste0("PC",kk), y=paste0("PC",jj), fill="ancestry"), shape=4, stroke=1)
		ggsave(paste0(OUTDIR,"/sample_PC",kk,"_PC",jj,".png"), g[jj,kk])
	}
	ggsave(paste0(OUTDIR,"/sample_PC1-4.png"), g)

}

PAIR<-read.csv(paste0(OUTDIR,"/relpairs.csv"), header=T, stringsAsFactors=F)
rownames(PAIR)<-with(PAIR,paste(ID1,ID2,sep=IDSEP))
PAIR$kinship <- 0.5*PAIR$kinship

# Identify Error Pairs
if (KINTOOL == 'king') {
	PAIR$DUPERR<-with(PAIR, TYPE=="Duplicate" & Kinship<MINDUP | TYPE!="Duplicate" & Kinship>=MINDUP)
	PAIR$POERR<-with(PAIR, TYPE=="Parent-Offspring"  & !(Kinship>=MIN1D & Kinship<MAX1D & IBS0<POIBS0))
	PAIR$SIBERR<-with(PAIR, TYPE=="Sib-Pair" & !(Kinship>=MIN1D & Kinship<MAX1D & IBS0>=POIBS0 ))
	PAIR$WFREL<-with(PAIR, TYPE=="Within-Family"  &  Kinship>=CRPREL & Kinship < MIN1D)
	PAIR$BFREL<-with(PAIR, TYPE=="Between-Family" &  Kinship>=CRPREL & Kinship < MIN1D)
	#PAIR$ISDUP<-with(PAIR, Kinship>=MINDUP)
	PAIR$WFPO<-with(PAIR, TYPE=="Within-Family"  &  Kinship>=MIN1D & Kinship<MAX1D & IBS0<POIBS0 )
	PAIR$WFSIB<-with(PAIR, TYPE=="Within-Family"  &  Kinship>=MIN1D & Kinship<MAX1D & IBS0>=POIBS0 )
	PAIR$BFPO<-with(PAIR, TYPE=="Between-Family" &  Kinship>=MIN1D & Kinship<MAX1D & IBS0<POIBS0 )
	PAIR$BFSIB<-with(PAIR, TYPE=="Between-Family" &  Kinship>=MIN1D & Kinship<MAX1D & IBS0>=POIBS0 )
	PAIR$FATAL<-with(PAIR, Kinship>=MAX1D & Kinship<MINDUP )
} else if (KINTOOL == 'peddy') {
	PAIR$DUPERR<-with(PAIR, TYPE=="Duplicate" & kinship<MINDUP | TYPE!="Duplicate" & Kinship>=MINDUP)
	PAIR$POERR<-with(PAIR, TYPE=="Parent-Offspring"  & !(kinship>=MIN1D & kinship<MAX1D & ibs0<POIBS0))
	PAIR$SIBERR<-with(PAIR, TYPE=="Sib-Pair" & !(kinship>=MIN1D & kinship<MAX1D & ibs0>=POIBS0 ))
	PAIR$WFREL<-with(PAIR, TYPE=="Within-Family"  &  kinship>=CRPREL & kinship < MIN1D  )
	PAIR$BFREL<-with(PAIR, TYPE=="Between-Family" &  kinship>=CRPREL & kinship < MIN1D )
	#PAIR$ISDUP<-with(PAIR, kinship>MINDUP)
	PAIR$WFPO<-with(PAIR, TYPE=="Within-Family"  &  kinship>=MIN1D & kinship<MAX1D & ibs0<POIBS0 )
	PAIR$WFSIB<-with(PAIR, TYPE=="Within-Family"  &  kinship>=MIN1D & kinship<MAX1D & ibs0>=POIBS0 )
	PAIR$BFPO<-with(PAIR, TYPE=="Between-Family" &  kinship>=MIN1D & kinship<MAX1D & ibs0<POIBS0 )
	PAIR$BFSIB<-with(PAIR, TYPE=="Between-Family" &  kinship>=MIN1D & kinship<MAX1D & ibs0>=POIBS0 )
	PAIR$FATAL<-with(PAIR, kinship>=MAX1D & kinship<MINDUP )
} else {
	stop(paste("Do not support ", TOOL), call.=FALSE)
}
ERRPAIR<-subset(PAIR,  DUPERR | POERR | SIBERR | WFREL | BFREL | WFPO | BFPO | WFSIB | BFSIB | FATAL) 


if(file.exists(paste0(OUTDIR,"/relpairs_err.csv")))
	file.remove(paste0(OUTDIR,"/relpairs_err.csv"))
if(nrow(ERRPAIR) > 0) {
	write.csv(ERRPAIR, paste0(OUTDIR,"/relpairs_err.csv"), quote=FALSE, row.names=FALSE)
}

if (NOPLOT == 0) {

	g<-ggplot(PAIR) + geom_point(aes(x=IBS0, y=Kinship, color=TYPE), shape="x", size=2, alpha=0.5) +
		labs(x="IBS0", y="Kinship Coefficient", color="Types of Pairs") 
	if (nrow(ERRPAIR) > 0) {
		g<-g+geom_point(data=ERRPAIR, aes(x=IBS0, y=Kinship, color=TYPE), size=2, shape=5) +
			#geom_text(data=ERRPAIR, aes(x=IBS0, y=Kinship), label=rownames(ERRPAIR), hjust=0, size=2, check_overlap=TRUE)
			geom_text_repel(data=ERRPAIR, aes(x=IBS0, y=Kinship), label=rownames(ERRPAIR), size=2)
	}
	ggsave.png(paste0(OUTDIR,"/relpairs_KinshipIBS0_king.png"), g, dpi=300)

	g<-ggplot(PAIR) + geom_point(aes(x=ibs0, y=kinship, color=TYPE), shape="x", size=2, alpha=0.5) +
		labs(x="IBS0", y="Kinship Coefficient", color="Types of Pairs")
	if (nrow(ERRPAIR) > 0) {
		g<-g+geom_point(data=ERRPAIR, aes(x=ibs0, y=kinship, color=TYPE), size=2, shape=5) +
		#geom_text(data=ERRPAIR, aes(x=ibs0, y=kinship), label=rownames(ERRPAIR), hjust=0, size=2, check_overlap=TRUE)
		geom_text_repel(data=ERRPAIR, aes(x=ibs0, y=kinship), label=rownames(ERRPAIR), size=2)
	}
	ggsave.png(paste0(OUTDIR,"/relpairs_KinshipIBS0_peddy.png"), g, dpi=300)

	g<-ggplot(PAIR) + geom_point(aes(x=Z1, y=Z2, color=TYPE), shape="x", size=2, alpha=0.5) +
		labs(x="Prob(IBD=1)", y="Prob(IBD=2)", color="Types of Pairs")
	if (nrow(ERRPAIR) > 0) {
		g<-g+geom_point(data=ERRPAIR, aes(x=Z1, y=Z2, color=TYPE), size=2, shape=5) +
			# geom_text(data=ERRPAIR, aes(x=Z1, y=Z2), label=rownames(ERRPAIR), hjust=0, size=2, check_overlap=TRUE)
			geom_text_repel(data=ERRPAIR, aes(x=Z1, y=Z2), label=rownames(ERRPAIR), size=2)
	}
	ggsave.png(paste0(OUTDIR,"/relpairs_IBD1_IBD2_plink.png"), g, dpi=300)

	g<-ggtern(PAIR) + geom_point(aes(x=Z0, y=Z1, z=Z2, color=TYPE), shape="x", size=2, alpha=0.5) + 
		xlab("P(IBD0)") + ylab("P(IBD1)") + zlab("P(IBD2)") + theme_showarrows() 
	if (nrow(ERRPAIR) > 0) {
		g<-g+geom_point(data=ERRPAIR, aes(x=Z0, y=Z1, z=Z2, color=TYPE), size=2, shape=5) +
		geom_text(data=ERRPAIR, aes(x=Z0, y=Z1, z=Z2), label=rownames(ERRPAIR), hjust=0, size=2, check_overlap=TRUE) +
		#geom_text_repel(data=ERRPAIR, aes(x=Z0, y=Z1, z=Z2), label=rownames(ERRPAIR))
		theme_nomask()
	}	
	ggsave(paste0(OUTDIR,"/relpairs_IBD_Ternary_plink.png"), g, dpi=300)

}