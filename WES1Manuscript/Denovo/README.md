### Discovery cohorts
#### Sample/cohort summary
We collected following information for four discovery cohorts (ASC, MSSNG, SSC, and SPARK), results are given previously in `202003_PubVars/AutismMergePatched/Samp_Summary.txt` and linked here as "TrioSampSummary.txt" (The table is displayed by the output of following command)

```bash
csvtk round -t -f 3,4,6 Discovery/TrioSampSummary.txt | csvtk csv2md
```

| Cohort   | Affected | M/F  | pMPX | nID   | pID  | Uncertain | Unaffected |
| :------- | :------- | :--- | :--- | :---- | :--- | :-------- | :--------- |
| ASC      | 4076     | 4.19 | 0.09 | 2705  | 0.43 | 0         | 351        |
| MSSNG    | 3363     | 3.96 | 0.35 | .     | 0.00 | 62        | 222        |
| SSC      | 2654     | 6.29 | 0.00 | 2562  | 0.41 | 0         | 2176       |
| SPARK    | 7015     | 3.96 | 0.25 | 6979  | 0.23 | 27        | 3017       |
| Combined | 16877    | 4.29 | 0.20 | 12245 | 0.32 | 89        | 5764       |

 	1. Affected: Unique number of affected and unaffected offspring. We used the same definion of pheno groups and sample inclusion/exclusion criteria in burden analysis. Note the dedup are only performed within sample group only, and we have confirmed with the DNV burden output (same for other sample groups).
 	2. M/F: Male-to-female ratio in cases. 
 	3. pMPX: Fraction of trios in families with known history (multiplex): for SPARK, we used FamType Simplex/Multiplex; for SSC, none have family history; for ASC, MPX is defined for trios who have other affected sibs in the same family based on FamID; for MSSNG, MPX is defined either from DB5 or more than one unquie cases in the same family. All trios without known family history are considered as simplex (no known family history).
 	4. pID: Fraction of low functioning cases (`ID` was used as proxy): for SPARK, we use self-reported cognitive impairment; for SSC, we use Vineland score<=70; for ASC, we used binary indicator of intellectual disability; for MSSNG, this information is not available. Only individuals who have non-missing information are used as denominator (nID is the number of individuals with ID info). 

Sample sizes are consistent with the output of the burden analysis  `202003_PubVars/AutismMergePatched/DNVBurden.sampsize.txt`.

We will also mark up analyzed trio samples from SPARK/SSC.

```bash
mkdir Discovery/TrioSamps
flatdb.pl --in ../../202003_PubVars/AutismMergePatched/SampTable.txt -f 'Pheno!="Uncertain"' --select IID,Pheno | table_intersect.pl -a ../Samples/SPARK/PedPhenoInfo.txt -a-select FamID,IID,Father,Mother -b - > Discovery/TrioSamps/SPARK_Trios.txt
awk 'NR>1{print $2, "x"; print $3, "x"; print $4, "x"}' Discovery/TrioSamps/SPARK_Trios.txt | sort -u > Discovery/TrioSamps/SPARK_SampList.txt
flatdb.pl --in ../../202003_PubVars/AutismMergePatched/SampTable.txt -f 'Pheno!="Uncertain"' --select IID,Pheno | table_intersect.pl -a ../Samples/SSC/PedPhenoInfo.txt -a-fsep "\t" -a-select FamID,IID,Father,Mother -b - > Discovery/TrioSamps/SSC_Trios.txt
awk 'NR>1{print $2, "x"; print $3, "x"; print $4, "x"}' Discovery/TrioSamps/SSC_Trios.txt | sort -u > Discovery/TrioSamps/SSC_SampList.txt
```

#### Burden

Burden results for four individual cohorts were collected from following sources (with column name remapping, merged to `Discovery/BurdenByCohort/Merged.txt` by Discovery/scripts/collect_burden_percohort.pl script):

* SPARK: `~/Projects/SPARK30K/201909_DNVs/Merged/DNVBurden`, 
  	 * Affected: Affected, Unaffected: UnaffNoID, Low func: AffwID, High func: AffwoID, No or unknown fam history: Simplex, Known fam history: Multiplex

* SSC: `~/Projects/SSC/DNVs/Burden/CodingDNVs`,  
  * Affected: Affected, Unaffected: Unaffected, Low func: AffVinelandLe70, High func: AffVinelandGt70
* MSSNG: `~/Projects/SPARK30K/202003_PubVars/MSSNGMerged/DNVBurden`
  * No or unknown fam history: Simplex2, known fam history: Multiplex
* ASC: `~/Projects/SPARK30K/202003_PubVars/ASCMerged/DNVBurden`
  * Affected: Case, Low func: CaseID, High func: CaseNoID

Cohort-specific burden will be collected and visualized to demostrate similar burdens across cohorts ("**Discovery/BurdenByCohort/**").  Fold enrichments across cohorts are similar if sequencing depth and family history are accounted for.

Burden for combined discovery cohort are previously calculated and can be found in `../202003_PubVars/AutismMergePatched`:

- Original result: `DNVBurden.burden.txt` is used to estimate PPV weights for DenovoWEST analysis
- Plot version: `DNVBurdenPlot.contra.txt` is the case-control burden comparison for plots (a related table `DNVBurdenPlotExclKnown.contra.txt` is used to for the same plot after excluding known genes)
- LoF version:  `DNVBurdenLoFReannoSummary.txt` summarize the influence of AF and pExt filters on de novo LoFs and Dmis.

The burden results from the above tables are summarized in `Discovery/BurdenMerged.xlsx`), we will generate following figures in `Discovery/BurdenMerged`:

1. We first summarize the burden of damaging de novo variants for combined cohorts across different gene sets.  We first considered genes sets defined by ExAC pLI>=0.5 or pLI<0.5. The plot shows the *de novo* burden and case-control rate difference in each gene set. We can found known ASD/NDD genes (the full list, 626 genes) explain 70% and 63% of LoF and Dmis rate differences among constrained genes, respectively. 
2. We also considered gene sets defined by gnomAD LOEUF deciles, and found de novo burden concentrates in the top 2 deciles. Known genes explain 71% and 59% LoF and Dmis rate difference in genes from top 2 deciles. 
3. We also compared the burden between affected with or without ID (ID is a proxy of lower functioning). And found affected with ID have higher de novo burden in constrained genes and especially in genes at top LOEUF decile. But this difference can be explained by known genes.
4. We also evaluated the influence of LoF prioritization method on de novo LoFs, we will later compare with inherited LoFs.

#### Gene set enrichment

We run dnenrich to test the enrichment of de novo variants in curated gene sets (the same set used in testing enrichment of inherited LoFs). The results are summarized in `../../201912_RareVars/TDT_Patch/GeneSetsDnEnrichSummary.txt`. We confirmed that most enriched gene sets can be explained by known genes. 

Gene set enrichment plot can be found in `Discovery/DnGeneSetEnrich.pptx`

#### Recurrent DNVs

We identified recurrent DNVs from ~17K autism case trios from discovery cohort after resolving sample duplicates. We will make use of previously generated list after applying patch for pExt and indicate known genes.

```bash
table_intersect.pl -a ../../202003_PubVars/AutismMergePatched/CodingDNVs_Recurrent_wMCRnPERs_pext.txt -a-fsep "\t" -b ../Genesets/Known/GeneTable.txt -b-select GeneID,Strict,Full -b-fsep "\t" -with-loj  > Discovery/RecurrentDNVs/AnnoPatched.txt
# Further check co-localization with standing LoFs
perl Discovery/scripts/colocalize_denovo_standingLoFs.pl
```

We perform manual curation of likely pathogenic DNVs among those recurrent in autism.

Recurrence between ASD and NDD can also be used to guide the interpretation of likely pathogenic DNVs. 

We identified 66 recurrent among 16877 unique autism cases among 15810 families, 16 of which overlap with other NDD DNVs. Among those autism recurrent DNVs, we identified 29 likely pathogenic based on our criteria, including 14 overlap with other NDD DNVs (suggesting they are pathogenic across conditions).  The properties of recurrent DNVs are compared with DDD (**Discovery/RecurrentDNVs.png**). The notable difference with NDD recurrent DNVs is: in autism, the fraction of pathogenic recurrent dnLoF is higher whereas  recurrent missense variants in genes with non-LoF mechanisms  is lower (P=0.03).


#### DenovoWEST/TADA
To identify gene significantly enriched by de novo variants, we applied both DenovoWEST and a modified version of TADA. We have summarize the gene level counts of different variant classes in `202003_PubVars/AutismMergePatched/DNVBurden_MergedGeneTab.txt`. PPV weights for DenovoWEST were calculated for LoF/Dmis/Bmis separately in LoF intolerant (pLI>=0.5) or other genes based on the burden results (`Discovery/DenovoWEST_weights.xlsx`). TADA priors were estimated using extTADA.

DenovoWEST and TADA results for the discovery sample are downloaded from server:

* `Discovery/DenovoWEST.txt` -- This is the full exome-wide scan in discovery cohort using DenovoWEST with autism-specific weights. The DenovoWEST-only results were reformatted to `Discovery/DenovoWEST.xlsx`. **10506** genes have DenovoWEST results, other genes do not have DNVs so p-value can be considered as 1.
* `Discovery/DenovoWEST/TopGenesHighSynOvsE.txt` -- Sensitivity analysis for top genes (P<1e-3) with SynOvsE>1.2, we noted that one gene (*HIST1H1E*) has unusual high SynOvsE, it will be removed from replication list.
* `Discovery/DenovoWEST/TopGenesSimplex.txt` -- Re-analysis of top genes (P<1e-3)  in families without known family history. This result will only be used as reference, focusing on simplex families does not necessarily increase power due to reduced sample size.	
* `Discovery/TADA_DNMisModel.txt` -- The output of modified TADA that stratify genes by pLI and also allow a fraction of disease genes to harbor missense variants only. A related file `TADA/DNMisModel_BFs.txt` contains BFs calculated for LoF and Dmis before transfomration (Note: This is from `Autism/TADA_DNMisModel2.txt` on Shannon).

The exome-wide scan results for the discovery cohort will be merged from gene level counts, DenovoWEST and TADA. We further filter the merged table by EntrezID because DenovoWEST results for genes with abnormal entrezIDs are enriched for artefacts. 

```bash
# Merge DenovoWEST and TADA for top genes
table_intersect.pl -a ../../202003_PubVars/AutismMergePatched/DNVBurden_MergedGeneTab.txt \
	-a-fsep "\t" -b Discovery/DenovoWEST.txt -b-select GeneID,AllObserved,AllExpected,pAllEnrich,MisObserved,MisExpected,pMisEnrich,MisEvents,MisDist,pMisCluster,pMisComb,pDenovoWEST -b-fsep "\t" -with-loj | \
	table_intersect.pl -a - -a-fsep "\t" -b Discovery/TADA_DNMisModel.txt -b-select GeneID,dnLGD,mutLGD,dnDmis,mutDmis,log10BF_All,PP_All,Qvalue \
	-b-fsep "\t" -with-loj | flatdb.pl --in - -fsep "\t" -f 'EntrezID!="." & EntrezID<10000000' |  table_intersect.pl -a - -a-fsep "\t" -b Discovery/TADA/DNMisModel_BFs.txt -b-select GeneID,log10BF_dnLGD,log10BF_dnDmis | csvtk sort -t -k PP_All:nr -k pDenovoWEST:n > Discovery/DenovoWEST_TADA/AllGenesMerged.txt
```

A total of 160 genes have pDenovoWEST<1e-3 (Top 100 of them are also ranked at top 100 by TADA) and will be subject to further test in replication cohort (below). We have further removed one gene (*HIST1H1E*) which have SynOvsE=3.8. Among top prioritized genes, 108 are known or possible NDD genes based on DDG2P.

The gene level info, per-gene de novo counts, DenovoWEST and TADA results were combined and reformatted in`Discovery/DenovoWEST_TADA.xlsx`. 

##### TADA

We modified original TADA model to incoporate gene level constraint and allow a fraction of disease genes to harbor only missense pathogenic variants. To support this approach, we compared goodness of fit of modified and original models, and found:

1. Modified model fits the observed de novo data better without introducing too many parameters (lower BIC)
2. The model captures the trend of fold enrichment and recurrence of de novo LoF and Dmis in constrained and non-constrained genes respectively

The TADA model fit and comparison between different models can be found in `Discovery/TADA` and summarize in `Discovery/TADA.xlsx` Note that we also applied modified TADA model to NDD data. The fitted model and calculated BFs will be used to compare LoF/Dmis enrichment among top genes between ASD and NDD.

We compared p-value and TADA's posterior probability (PP), and found top 100 TADA genes were also have P<1e-3 by DenovoWEST. Results are merged with DenovoWEST results in `Replication/DenovoWEST_TADA/AllGenesMerged.txt` and graphically summarized in `Discovery/DenovoWEST_TADA`.

#### Likely pathogenic DNVs

To calibrate filters for inherited LoFs, we perform comprehensive annotations in selected genes and compare de novo vs inherited (In inherited part).


### Replication
#### Sample/cohort summary

We included trios from SPARK WES2, WES3, and WGS for replication. We used released data for WES2 (2020-06), WGS1 (2020-02), and WGS2-3 (2020-12), and unreleased data for WES3 (2020-10 from Regeneron). Sample phenotype and family type were collected from master table for released data and from phenotype database V4 for unreleased data (or V3 for family type). 

Note that de novo from AGRE were also available, but we found it have substantial overlap with MSSNG but difficult to find out the overlapping sample due to changes in sample IDs. And also because AGRE families are mostly multiplex, we did not include AGRE as part of the replication cohort. 

To summarize sample info, we will add cohort level info, and will further collect ID (cognitive impairment or DevID) info from DB V4. To determine total number of unique samples, we further keep one of twin pairs for each sample group and remove likely overlapping samples with discovery cohort. Samp level info for the replication cohort can be found in `Replication/SampInfo`. 

To tally cohort level summary, we first run `Replication/scripts/append_sampinfo.pl` to create sample table in the same format as discovery cohort, then run `Replication/scripts/sum_rep_cohorts.pl` generate basic cohort level summary (results are appended to `DnTrioSummary.xlsx`) and shown below:  

| Cohort     | Affected | M/F  | pMPX | nID   | pID  | Uncertain | Unaffected |
| :--------- | :------- | :--- | :--- | :---- | :--- | :-------- | :--------- |
| Discov     | 16877    | 4.29 | 0.20 | 12245 | 0.32 | 89        | 5764       |
| RepWES     | 4137     | 3.62 | 0.21 | 4104  | 0.24 | 16        | 831        |
| RepWGS     | 2088     | 4.33 | 0.14 | 2081  | 0.23 | 12        | 1663       |
| Rep        | 6188     | 3.85 | 0.18 | 6148  | 0.24 | 28        | 2472       |
| CombAutism | 23053    | 4.16 | 0.19 | 18392 | 0.29 | 117       | 8236       |
| CombNDD    | 54588    | 1.98 |      |       |      |           |            |

#### Burden and DenovoWEST/TADA

Gene level burden have been summarized after accounting for duplicated sample. We perform DenovoWEST meta-analysis first using weights derived from discovery cohort, we also run TADA on top genes using priors from discovery cohort.

The following files are downloaded from server:

* `Replication/Burden/Meta_Autism.genetab.txt` -- Gene level summary of top genes in discovery+replication autism cohorts. (Meta_Autism2.genetab.txt contains results for prioritized top TD genes.)
* `Replication/Burden/Meta_NDD.genetab.txt` -- Gene level summary of top genes in full autism cohorts plus NDD cohorts (mostly from DDD cohort).
* `Replication/DenovoWEST/Meta_Autism.txt` -- DenovoWEST meta-analysis of top genes in discovery+replication  autism cohorts. (Note: `Meta_Autism2.txt` is the results for prioritized top TD genes).
* `Replication/TADA/DNMisModel.txt` -- TADA meta-analysis of top genes in discovery+replication cohorts. (DNMisModel2.txt contains the result for prioritized TD genes).
* `Replication/TADA/DNMisModel_BFs` -- The same as above with BFs calculated for Dmis and LoF.
* `Replication/DenovoWEST/Meta_AutismSimplex.txt` --  Meta-analysis of all autism cases who have no known family history. (`Meta_AutismSimplex2.txt` contains the results for prioritized top TD genes).
* `Replication/DenovoWEST/SPARK+SSC_Autism.txt` -- DenovoWEST analysis of top genes in SSC+SPARK cohorts only. This is aimed to remove the effect of sample ascertainment in ASC or MSSNG.
* `Replication/DenovoWEST/SPARK+SSC_AutismSimplex.txt`-- Same as above, but only include simplex cases.
* `Replication/DenovoWEST/RepOnly_Autism.txt` -- DenovoWEST analysis of top genes in replication cohorts only.
* `Replication/DenovoWEST/RepOnly_AutismSimplex.txt` -- DenovoWEST analysis of top genes in autism cases without known family history from replication cohorts only.
* `Replication/DenovoWEST/Meta_NDD.txt` -- DenovoWEST meta-analysis of top genes in combined ASD-NDD cohort.
* `Replication/DenovoWEST/Meta_NDDSimplex.txt `--  Same as above, but excluding autism cases with known family history.

Then merge results for top genes to `Replication/DenovoWEST_TADA/TopGenesMerged.txt` (by Replication/scripts/merge_topgene_tabs.pl). We have reformatted the results in `Replication/DenovoWEST_TADA_TopGenes.xlsx`

Note: we also collected DenovoWEST/TADA results for prioritized TD genes, they will appear in the second sheet of the above file.



Summary of the results:

1. A total of 46 genes are exome-wide significant in discovery cohort (45 are also EWS in discovery+replication, the remaining 1 is at borderline), 61 are exome-wide significant in discovery+replication. Most EWS genes are known NDD or autism genes, RFX3 and RORB are not included in DDG2P (2020-02), and have converging evidence from recent literature. MARK2 is the novel genes not reported before.
2. Using TADA PP>0.95 prioritize additional genes that do not reach EWS, including novel candidates TAF4, FUBP3, SRPRA, TNPO3, PRKAR1B and PAPOLG. <- Notet that only PRKAR1B have strong replication evidence.
3. Most EWS genes (59 of 61) in autism discovery+replication also have evidence with P<0.01 in DDD (mainly NDD), except SMARCC2 and BRSK2 which is likely due to sample ascertainment or sequencing depth. But SMARCC2 was later found dubious in mega-CC analysis. Most were also significant in SPARK+SSC subset, except for KCNQ3, which is likely due to difference in sample ascertainment.
4. Meta-analysis of ASD (discovery+replication) and NDD reveals a total 100 EWS genes, including all 62 EWS genes in ASD analysis.  Two novel genes prioritized by TADA in autism analysis now surpass EWS: PRKAR1B and TNPO3, SRPRA is at the borderline. Note that TNPO3 does not have significant higher LoF rate than population controls, but de novo signals are not EWS for missense variants only.
5. SMAD6 is EWS in autism+NDD meta-analysis and not known ASD/NDD gene. But this gene have high LoF OvsE=2, so was removed.
6. The relative contribution of LoF and Dmis can be quantified by BF[LoF]/BF[Dmis]. We compared this quantitify for 82 EWS genes that have PP>0.5 in both autism and NDD. A few genes with signficant different BF[LoF]/BF[Dmis] in autism and NDD can be noted.

We will seek further replication evidence for the above genes in de novo + TDT + pseudo CC meta-analysis.

```bash
# Extract gene lists for future analysis
# Top 61 genes that are exome-wide signif from ASD de novo meta-analysis
# 58 were in constrained genes list
flatdb.pl --in  Replication/DenovoWEST_TADA/TopGenesMerged.txt -fsep "\t" -f 'pDenovoWEST_Meta!="." & pDenovoWEST_Meta<0.0000013' --select GeneID,CytoBand,ExACpLI,LOEUFbin,HGNC > Replication/DenovoWEST_TADA/TopGenes_ASDMeta_ExomeSignif.txt
# Top 100 (98, excl TNPO3/SMAD6) genes that are exome-wide signif
flatdb.pl --in  Replication/DenovoWEST_TADA/TopGenesMerged.txt -fsep "\t" -f 'pDenovoWEST_ASDNDD!="." & pDenovoWEST_ASDNDD<0.0000013 & HGNC!~"^(SMAD6|TNPO3)"' --select GeneID,CytoBand,ExACpLI,LOEUFbin,HGNC > Replication/DenovoWEST_TADA/TopGenes_ASDNDDMeta_ExomeSignif.txt
# All top genes carried forward for replication
flatdb.pl --in Replication/DenovoWEST_TADA/TopGenesMerged.txt -fsep "\t" -f 'HGNC!="HIST1H1E"' --select GeneID,CytoBand,ExACpLI,LOEUFbin,HGNC > Replication/DenovoWEST_TADA/TopGenes_ForReplication.txt
```

Some genes have de novo LoFs in controls. They will be extracted for manual review

```BASH
table_intersect.pl -a ../Inherited/Meta/CandGenes.txt -a-select Known,Ascertainment,GeneID,CytoBand -a-fsep "\t" -b Replication/Burden/Meta_CaseCtrl.genetab.txt -b-fsep "\t" | flatdb.pl --in - -fsep "\t" -f 'MaleControl_LoF_Count>0 | FemaleControl_LoF_Count>0' > Replication/BurdenCaseCtrl.txt
# For comparison, also extract from discovery cohort
table_intersect.pl -a Replication/BurdenCaseCtrl.txt -a-select GeneID -a-fsep "\t" -b ../../202003_PubVars/AutismMergePatched/DNVBurdenCaseCtrl.genetab.txt -b-fsep "\t" -with-loj > Discovery/BurdenCaseCtrl.txt
```

#### Final data set

Reformat  the final data set used in publication.

Note we have downloaded DNVs from replication cohort to `Replication/DNVs`.

```bash
# Reformat ASD samp tabs => Data/ASD_{Discov,Rep}_Trios.txt
perl Data/scripts/refmt_Samptabs.pl
# Reformat ASD/NDD DNVs => Data/ASD_{Discov,Rep}DNVs.txt, Data/NDD_DNVs.txt
perl Data/scripts/refmt_DNVtabs.pl
# NDD Trios
cp ../../202003_PubVars/NDD/SampTable_wPheno.txt Data/NDD_Trios.txt 

```

