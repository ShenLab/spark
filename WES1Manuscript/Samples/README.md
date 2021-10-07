### SPARK
SPARK 30K sample was formally published for the first time. So we need to give more detailed information about detailed sample level info that were used in the analysis.
#### Basic Pedigree/Phenotype Information
We first extract key phenotype data used in other analysis and store a simplified table for summary: `SPARK_/PedPhenoInfo.txt` (by **SPARK/scripts/extract_SPARK_pheno.pl**). This table include **Pilot+30K** samples, regardless of data structure. All phenotype info are based on 2019_Version3 (may have small diff with 2020_Version4). Sample included in different types of analyses will be listed and summarized in analysis specific directories, they can be merged with the sample table at later stage. The following columns are in the sample table

 * ASD is based on ASD flag, presumed to be derived mostly from clinical diagnosis.
 * Diagnosis: one of ASD, Asperger, Autism, or PDD(NOS), available for >99.5% cases
 * ID is calculated from CogImpair at enrollment (for ASD cases only) and DevID from BMS questionaire (available for 85% total participates). 
   * Calculated ID info is available (not unknown) for >99% of ASD offspring, ~90% of nonASD offspring, and ~85% of parents.
 * Full scale IQ and Vineland (v3) score is available for a subset of individuals (pulled from latest v5 database). Note: IQ is only available for affected. 
   * IQ and Vineland score will be used to demostrate self-reported cognitive impairment is a proxy of adpative behavior.
 * FamID is the family IDs and connected by 1st deg cryptic relpairs if avaialble, FamIDOrig is the original ID (this is used to define relation w.r.t to proband). 
 * Family type is based on the presence or absence 1st affected relative pairs (no matter enrolled or not) in the family, this is defined on each family (based on FamID) not on individual.
 * Sex is genetically verified sex.
	* Role is the **relation with proband** in original family (based on FamIDOrig): one of Father, Mother, Proband, Sibling, Maternal or Paternal Half_Sibling. Pairwise relationships have been verified genetically. 
	 * Note: half sib-ships are provided in the database and have been confirmed to be  consistent with kinship estimates, but we did not further verify by mitocondria inheritance, because not all parents are available. But the meaning of Role can be confusing if proband is not included in sequencing, so in the future we may not use this column.
	 * Note2: If proband is sequenced and the sib is genetically verified identical to proband, we will modified role to MZTwin_Sibling to indicate the twin pair. 

##### Twin Pairs

MZ twin pairs should also be extracted: `SPARK/MZTwins.txt` (by **SPARK/scripts/extract_MZtwins.pl**). For those exist in the final table, **64 MZ twin pairs**. 53 are concordant in ASD phenotype (51 affected, 2 unaffected).

* Note: Twin pairs will not appear in Role if neither of them is proband.

DZ twins can be identified from full sib pairs with the same DoB, but information may not be complete, so we will not reformat them in sample list. 

##### Summary

Basic phenotype from above table are summarized by in: `SPARK/PedPhenoInfo/Summary.txt` (by **SPARK/scripts/sum_SPARK_pheno.pl**).

A few observations:

* Comparing affected between offspring and parents, parents are less likely have ID.

* Comparing sex ratio between simplex vs multiplex, we can see reduced M/F ratio in affected offspring from multiplex families (~3:1) and in affected parents (<1).
  * Note: more affected mothers than fathers in affected parents is likely due to ascertainment bias

* There are also a small number of non-ASD sample have ID
  * Note: we will define "controls" in analysis as those have no ASD and no ID).

#### Parent education
We used parent education level in calibrating PRS performance and in testing the effect of rare LGD variants. So we should make a summary of Edu distribution for samples included in the analysis: note that we mapped reported education level to *ISCED 1997 code* as follows: 

* graduate_or_professional_degree => 6, 
* baccalaureate_degree => 5,  
* current_college, associate_degree, trade_school, some_college => 4, 
* high_school_graduate, ged_diploma => 3, 			
* some_high_school => 2, 
* did_not_attend_high_school => 1. 

We have incoporated this to the [sample info table](####Basic Pedigree/Phenotype Information). In PRS analysis, only samples from 30K release with SNP genotype are used.

Basic Pedigree/Phenotype Information

#### QC Plots (30K sample only)

##### Pairwise relationships

Verifying pedigree relationships using genome-wide genotypes for all released sample.
Plotting pedigree relationship vs different genotype relatedness metrics for samples with array genotypes (`SPARK_SampSum/Relpairs_SNPGeno.png`) and using exome genotypes (`SPARK_SampSum/Relpairs_ExomeGeno.png`). Some UPD samples will appear as outlier PO pairs in the plot.

Cryptic related families was summarized in `SPARK_Merged.pptx`. Family ID will be merged for those merged families.

##### Sex check
Verify sex based on chrX/Y LRR or read depth signals.
One figure will show scattterplot for SNP array and exome: `SPARK_SampSum/Sexcheck.png`, will higlight sex chromsome abnormalities outliers manually later. 

	table_intersect.pl -a ../../201906_SampQC/share/20190729_ArrayQC/chrom_meanLRR.txt -a-select IID,Std:LRRSD,chrX,chrY -b SPARK_SampTab/PedPhenoInfo.txt -b-select IID,Sex > SPARK_SampSum/Sexcheck_SNParray.txt
	table_intersect.pl -a ../../201908_BamQC/BamStats/201908/chrom_depth.csv -a-select IID,NDP_chrX,NDP_chrY -b SPARK_SampTab/PedPhenoInfo.txt -b-select IID,Sex > SPARK_SampSum/Sexcheck_Exome.txt

##### Ancestry prediction 
Two figure for PCA projection to K1G+HGDP ancestral populations colored by predicted ancestral components: `SPARK_SampSum/PCA_SNPGeno.png` and `SPARK_SampSum/PCA_ExomeGeno.png`.
We used Prob(EUR) cuoff of 0.85 to predict European populations.

	flatdb.pl --in ../../201906_SampQC/PopPCA/SPARK_K1G+HGDPNAT/out/SPARK30K.Ancestry.txt -f 'Label!="SPARK30K"' > SPARK_SampSum/PCA_SNPGeno.txt
	table_intersect.pl -a ../../201906_SampQC/PopPCA/SPARK_K1G+HGDPNAT/out/SPARK30K.Ancestry.txt -b SPARK_SampTab/PedPhenoInfo.txt -b-select IID | awk 'NR>1' >> SPARK_SampSum/PCA_SNPGeno.txt
	flatdb.pl --in ../../201906_SampQC/PopPCA/SPARK_K1G+HGDPNAT_Exome/out/SPARK30K.Ancestry.txt -f 'Label!="SPARK30K"' > SPARK_SampSum/PCA_ExomeGeno.txt
	table_intersect.pl -a ../../201906_SampQC/PopPCA/SPARK_K1G+HGDPNAT_Exome/out/SPARK30K.Ancestry.txt -b SPARK_SampTab/PedPhenoInfo.txt -b-select IID | awk 'NR>1' >> SPARK_SampSum/PCA_ExomeGeno.txt
	Rscript SPARK_SampSum/scripts/PCA.R	

### SSC
SSC samples have mostly been published before. Some published samples do not have inhouse sequence level data. We generated the sample table to include all published and reprocessed samples that we included in analysis.
#### Basic pedigree / phenotype information
The simplified sample table `SSC/PedPhenoInfo.txt` was prepared by (**SSC/scripts/extract_SSC_pheno.pl**). We have following fields in the table:

* ASD is the autism affection status, True for proband, False for others.
* Diagnosis: CpeaDx clinical diagnosis one of 1-autism, 2-asd, 4-asperger, available for >96% cases	
* FullScaleIQ, VinelandIICompStdScore: Full scale IQ and full-scale Vineland (v2) ABC scores.
* FamID is the first five digits extracted from IID.
* All families are "presumed" simplex, although some of them are labeled as multiplex later, so we do not have FamType for SSC.
* Sex is genetically verified sex, Role is the relation with proband in the family.
* WES,WGS are indicator if WES or WGS data passed QC and processed in-house.

Note that only one MZ twin pair was included with ID suffix p2 in the sample table.

Need to summarize the platform for SSC sample (we will do this manually).

Also append columns to indicate if a sample is included certain type of analysis.

### Other Cohorts

We also included MSSNG (db6) and ASC (2014+2018) trios in the discovery cohort for de novo analysis. Variants were collected from public sources.

For replication of de novo vatiants, SPARK replication cohorts will be summarized in de novo directory.

For replication of inherited variants, we make use of AGRE, with variants contributed by collaborators. 

#### MSSNG

We used selected samples from DB6 release (which is a superset of Yuen et al. 2017 DB4) for de novo variants. Only sample with confirmed autism/unaffected phenotypes, used by Yuen et al. 2017 for de novo calling, or have less than 250 total DNV calls are included:  ``` '(DNVCount == "." | DNVCount < 250 | DB4DNV=="x") & Pheno!="0" & Pheno!="."'.``` MSSNG trios are listed in `Other/MSSNG_TriosInfo.txt`

In additional to sample availability, the following fields are collected and included in the sample table:
 * FamID is the designated family ID from latest DB6.
 * FamType: family type (simplex or multiplex) is synthesized from Yuen 2017 (DB4), DB5, or observed number of cases in the same family (available for ~65% of all families).
 * Sex is based on the latest DB6 annotation (we did not check for consistency across all publications)
 * DNASource,Platform: DNA source and sequencing platforms for each sample, will be useful to an initial evaluatation of the variant call quality.
 * ASD is the autism affection status, we only include sample that we can collect phenotype from literature or infer as proband based on sample IDs (indicated by PhenoSrc). Note: sample with pheno==0 have uncertain or non-autism phenotypes and are not included.

#### ASC

We merged published samples in ASC2014 and ASC2018 studies: (`../../202003_PubVars/ASCMerged/SampTable.txt`). Identified 6 pairs of likely dup/twins (`../../202003_PubVars/ASCMerged/Samp_LikelyDupPairs.txt`).

#### AGRE (not included)

We will reformat the table files contributed by collaborators: 3926 samples from 859 families (`Other/AGRE_SampInfo.txt`) including 79 MZ twin pairs (`Other/AGRE_MZTwins.txt`).



​	

