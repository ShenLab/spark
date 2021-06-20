# Calculate BF and PP using TADA model
# 1. Input variant count table: GeneID,dnCount,mutRate,CaseCount.NCase,CtrlCount.NCtrl,...
# 	After first column (gene ID), every two columns represent counts and rate (or control count) for a variant class
#   Note: For case/control data, we will pack sample size to both column names (".SampSize")
#		  For de novo data, we pack sample size to the count column.
# 2. Table for prior parameters, (GeneID),Pi,RRMean1,Beta2,RRMean2,Beta2
#	After first column (fraction of disease gene), every two columns represent RR/Beta for a variant class 
#   in the same order as appearing in the input variant count table.
#   It is possible to have multiple rows in the table, e.g. when parameters are sampled from MCMC.
#	In such case, BF/PP calculation will be the average (numerical integration) over all samples of the parameters.
#   We can also specify gene-specific prior parameters by adding GeneID to the first column
#   in such cases, we require the priors should exist for each gene in the count table.
# 3. We currently does not allow sample size to vary across genes, but it can be achieve by proper scaling
#   of mutation rates for de novo data. In the future, we may allow variable sample sizes by further 
#   including a sample size table
# 4. The output is Bayes factor for each variant class and combined BF for all variant classes.
# Update 202005:
# The tradition TADA integrate evidence from multiple types of variants for the same set of disease genes. But we know 
# there are disease genes targeted by missense variants only. To accomondate this, we modified TADA framework to model
# the disease genes as 2 or more classes. Each additional class represent te genes targeted by certain groups of variants
# only. Each disease gene class will be associated with a prior prob of being in this class, and prior RR/beta for 
# the effect size. When a set of gene was targeted by missense variants only, the prior for LGD variants should be set
# to 1 instead of missing. So we now allow multiple set of priors from input.
# Update 202006: we have changed natural log to log10. Log10_BF has much more interpretable than log_BF.


# Bayes factor for de novo variants: account for NA
bayes.factor.denovo <- function(obs, N, mu, gamma.mean, beta) {
	marg.lik0 <- dpois(obs, 2*N*mu)
	marg.lik1 <- dnbinom(obs, gamma.mean*beta, beta/(beta+2*N*mu))
	BF <- marg.lik1/marg.lik0
	BF[is.na(BF)]<-1
 	return(BF)
}

# Bayes factor for case-control data (adapted from extTADA)
# This is the vectorized version, account for NA.
bayes.factor.cc <- function(x.case, x.control, Nsample, gamma.meanCC, betaCC, rhoCC, nuCC) {
	if(is.na(gamma.meanCC) || is.na(betaCC)) {
		return(rep(1, length(x.case)))
	}

    gAll <- range(rgamma(10000, gamma.meanCC*betaCC, rate = betaCC))
    gLower = gAll[1]; gUpper = gAll[2]

    altCC <- apply(cbind(x.case, x.control), 1, function(y){
        x2 <- list(ca = y[1], cn = y[2])
        if(is.na(y[1]) || is.na(y[2])) {
        	return(1)
        }
        evidence.alt.cc3 <- function(x = x2, N = Nsample, gamma.mean = gamma.meanCC, beta = betaCC,
                                     rho1 = rhoCC, nu1 = nuCC) {
            bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
            fCase <- function(gGamma) {
                dnbinom(x$ca, size = rho1 + x$cn, prob = (N$cn + nu1)/(N$cn + nu1 + N$ca*gGamma))*dgamma(gGamma, gamma.mean*betaCC, rate = betaCC)
            }
            bCase <- log(integrate(fCase, lower = gLower, upper = gUpper, stop.on.error = FALSE)$value)
            return(exp(bCase+bControl))
        }
        t1 <- evidence.alt.cc3()
        return(t1)
    })

    nullCC <- apply(cbind(x.case, x.control), 1, function(y, rho1 = rhoCC, nu1 = nuCC, N = Nsample){
        x <- list(ca = y[1], cn = y[2])
         if(is.na(y[1]) || is.na(y[2])) {
        	return(1)
        }
        bControl <- log(dnbinom(x$cn, size=rho1, prob = nu1/(N$cn + nu1)))
        bCase <- log(dnbinom(x$ca, size=rho1+x$cn, prob=(N$cn+nu1)/(N$cn+nu1+N$ca)))
        t1 <- exp(bCase+bControl)
        return(t1)
    })

    tempBF <- altCC/nullCC #ifelse((x.case == 0) & (x.control == 0), 1, altCC/nullCC)
    return(tempBF)
}


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
	stop("Must provide input variant count table, one or more prior tables, and output file name!");
}

# Reading variant count table
VarCount<-read.table(args[1], header=TRUE, na.string=".", as.is=TRUE)
if(ncol(VarCount) %% 2 != 1) {
	stop("Incorrect number of columns in variant count table!")
}

# Read one or more prior sets
PriorSets<-list()
for(ii in seq(2, length(args)-1)) {
	Priors<-read.table(args[ii], header=TRUE, na.string=".", as.is=TRUE)
	if (ncol(VarCount) != ncol(Priors) && ncol(VarCount) != ncol(Priors)-1) {
		stop(paste0("Incorrect number of columns in priors file: ", args[ii]))
	}
	PriorSets[[ii-1]]<-Priors
}
# When multiple set of priors are available, we need to check they have the same number of rows
NPriorRows<-nrow(PriorSets[[1]])
if(length(PriorSets) > 1) {
	for(ii in seq(2, length(PriorSets))) {
		if(nrow(PriorSets[[ii]]) != NPriorRows) {
			stop("Inconsistent number of rows among prior tables!")
		}
	}
}

OutFile<-args[length(args)]

# Determine the number of classes and their types
VarType=strsplit(names(VarCount), '.', fixed=TRUE)
OutFields<-vector()
for(ii in seq(1, (length(VarType)-1)/2)) {
	#OutFields<-c(OutFields, paste0("logBF_", VarType[[2*ii]][1]))
	OutFields<-c(OutFields, paste0("log10BF_", VarType[[2*ii]][1]))
}

if (ncol(VarCount) == ncol(Priors)) {
	# When columns in VarCount are the same as Prior table
	# Go through each variant class and aggregate evidence
	# 
	AvgLogBF<-rep(0, nrow(VarCount))
	AvgPP<-rep(0, nrow(VarCount))
	AvgLogBFvar<-matrix(0, nrow=nrow(VarCount), ncol=(ncol(VarCount)-1)/2)

	# When there are multiple rows in the prior table, they will be taken as multiple samplings
	# from MCMC
	for(jj in seq(1, NPriorRows)) {
		# Go-through each prior set
		BFall<-list()
		BFvar<-list()
		# pi0 is now summed over all prior sets
		pi0<-sum(unlist(lapply(PriorSets, function(a) a[jj, 1] )))
		if (pi0>1) {
			stop("Incorrect prior prob of disease genes!")
		}
		for(kk in seq(1, length(PriorSets))) {
			Priors<-PriorSets[[kk]]
			BFall[[kk]]<-rep(1, nrow(VarCount))
			BFvar[[kk]]<-matrix(0, nrow=nrow(VarCount), ncol=(ncol(VarCount)-1)/2)

			# Within each prior set
			# Go-through each variant type to aggregate evidence in terms of BF
			for(ii in seq(1, (ncol(VarCount)-1)/2)) {
				# Denovo 
				#message(paste("Processing", VarType[[2*ii]][1]))
				if (length(VarType[[2*ii]]) == 2 && length(VarType[[2*ii+1]]) == 1) {
					Ntrio<-as.numeric(VarType[[2*ii]][2])
					bfVar<-bayes.factor.denovo(obs=VarCount[,2*ii], N=Ntrio, mu=VarCount[,2*ii+1],
											gamma.mean=Priors[jj,2*ii], beta=Priors[jj,2*ii+1])

				} else if (length(VarType[[2*ii]]) == 2 && length(VarType[[2*ii+1]]) == 2) {
				# case-control
					Ncase<-as.numeric(VarType[[2*ii]][2])
					Nctrl<-as.numeric(VarType[[2*ii+1]][2])
					Nu<-200
					Rho<-Nu*mean(VarCount[,2*ii]+VarCount[,2*ii+1], na.rm=TRUE)/(Ncase+Nctrl)
					#cat(Ncase, Nctrl, Nu, Rho, fill=TRUE)
					bfVar<-bayes.factor.cc(Nsample=list(ca=Ncase, cn=Nctrl), 
										x.case=VarCount[,2*ii], x.control=VarCount[,2*ii+1],
										gamma.meanCC=Priors[jj,2*ii], betaCC=Priors[jj,2*ii+1],
										rhoCC=Rho, nuCC=Nu)
				} else {
					stop(paste0("Incorrect type of variant ", VarType[[2*ii]][1]))
				}
				BFall[[kk]]<-BFall[[kk]]*bfVar
				BFvar[[kk]][,ii]<-bfVar
			}
			BFall[[kk]]<-BFall[[kk]]*Priors[jj,1]/pi0
			BFvar[[kk]]<-BFvar[[kk]]*Priors[jj,1]/pi0
		}
		# Now convert to logBF 
		#logBFall<-log(Reduce('+', BFall))
		#logBFvar<-log(Reduce('+', BFvar))
		logBFall<-log10(Reduce('+', BFall))
		logBFvar<-log10(Reduce('+', BFvar))

		# Convert BF to PP (require Pi0)
		#PP<-pi0*exp(logBFall)/((1-pi0)+pi0*exp(logBFall))
		PP<-pi0*10^(logBFall)/((1-pi0)+pi0*10^(logBFall))

		# Take average over different samplings of parameters
		AvgLogBF<-AvgLogBF+(logBFall-AvgLogBF)/jj
		AvgPP<-AvgPP+(PP-AvgPP)/jj
		AvgLogBFvar<-AvgLogBFvar+(logBFvar-AvgLogBFvar)/jj

		# Write out BF for each variant class and BF/PP for all classes
		BFPP<-as.data.frame(cbind(VarCount[,1],AvgLogBFvar,AvgLogBF,AvgPP))
	}

} else {
	# Use gene-specific priors, require that all genes must have prior specified
	# gene-specifc prior is aimed to replicate ASC's results
	# There is no way to incoporate uncertainties
	# We also do not support multiple priors
	if(length(PriorSets) > 1) {
		stop("We do not support multiple gene-specific prior sets!")
	}
	rownames(Priors)<-Priors[,1]
	# Go-through each gene to find gene specific priors
	for(jj in seq(1, nrow(VarCount))) {
		GeneID<-VarCount[jj,1]
		if(! GeneID %in% rownames(Priors)) {
			stop("Cannot find priors for gene ", GeneID)
		}
	}
	GenePriors<-Priors[VarCount[,1], -1]

	logBFall<-rep(0, nrow(VarCount))
	logBFvar<-matrix(0, nrow=nrow(VarCount), ncol=(ncol(VarCount)-1)/2)
	
	# Go-through each vriant class
	for(ii in seq(1, (ncol(VarCount)-1)/2)) {
		if (length(VarType[[2*ii]]) == 2 && length(VarType[[2*ii+1]]) == 1) {
			Ntrio<-as.numeric(VarType[[2*ii]][2])
			bfVar<-bayes.factor.denovo(obs=VarCount[,2*ii], N=Ntrio, mu=VarCount[,2*ii+1],
									gamma.mean=GenePriors[,2*ii], beta=GenePriors[,2*ii+1])

		} else if (length(VarType[[2*ii]]) == 2 && length(VarType[[2*ii+1]]) == 2) {
		# case-control
			Ncase<-as.numeric(VarType[[2*ii]][2])
			Nctrl<-as.numeric(VarType[[2*ii+1]][2])
			Nu<-200
			Rho<-Nu*mean(VarCount[,2*ii]+VarCount[,2*ii+1])/(Ncase+Nctrl)
			# If priors are gene-specific, we cannot use vectorized version
			bfVar<-rep(0, nrow(VarCount))
			for(jj in seq(1, nrow(GenePriors))) {
				if(rownames(GenePriors)[jj] != VarCount[jj,1]) {
					stop("Gene prior table is not aligned with variant count table: ", jj)
				}
				bfVar[jj]<-bayes.factor.cc (Nsample=list(ca=Ncase, cn=Nctrl), 
									x.case=VarCount[jj,2*ii], x.control=VarCount[jj,2*ii+1],
									gamma.meanCC=GenePriors[jj,2*ii], betaCC=GenePriors[jj,2*ii+1],
									rhoCC=Rho, nuCC=Nu)
			}
			
		} else {
			stop(paste0("Incorrect type of variant ", VarType[[2*ii]][1]))
		}
		#logBFall<-logBFall+log(bfVar)
		#logBFvar[,ii]<-log(bfVar)
		logBFall<-logBFall+log10(bfVar)
		logBFvar[,ii]<-log10(bfVar)

		# Convert to PP
		#PP<-GenePriors[,1]*exp(logBFall)/(1-GenePriors[,1]+GenePriors[,1]*exp(logBFall))
		PP<-GenePriors[,1]*10^(logBFall)/(1-GenePriors[,1]+GenePriors[,1]*10^(logBFall))
		BFPP<-as.data.frame(cbind(VarCount[,1],logBFvar,logBFall,PP))
	}
}


if(ncol(BFPP) != length(OutFields)+3) {
	stop("Incorrect number of columns for final BFPP table!")
}
#names(BFPP)<-c("GeneID", OutFields, "logBF_All", "PP_All")
names(BFPP)<-c("GeneID", OutFields, "log10BF_All", "PP_All")
write.table(BFPP, file=OutFile, quote=FALSE, sep="\t", na=".", row.names = FALSE, col.names = TRUE)

