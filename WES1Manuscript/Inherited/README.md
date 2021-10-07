### Discovery cohorts

#### Sample/cohort summary

Sample size for TDT, pseudo case-control analysis were summarized in the following spreadsheet. 

`../../201912_RareVars/TDT_Patch/SampSizes.xlsx`

#### LoF filters

Before burden analysis we first need to prioritize inherited LoFs. This is done by evaluating pExt/lofee filters on de novo and inherited variants. For autosomal LoFs, applying pExt>0.1 and URLowAF filter, 45% standing LoFs in constrained genes were removed that do not contribute to TD signal. In comparison, less than 5% of de novo variants were removed.  Results were summarized in `TDBurden/pExtEval.pptx`

Note: that the above pExt threshold is tuned to maximize sensitivity, for individuals the optimal threshold is likely different. For de novo LoF enriched genes, we can contrasting the pattern of de novo and standing LoFs. 

In addition, we also found from the initial QC passed variant calls, up to 5% annotated LoF calls were part of nonLoF MNVs or frame-restoring indels in constrained genes. Results were summarized in `TDBurden/NonLoFsSummary.xlsx`

pExt filter will be re-evaluated later on known de novo LoF enriched genes. 

#### Burden analysis

Burden of discovery cohort were collected from the following files:

- `../../201912_RareVars/TDT_Patch/DiscovLoFTDBurden_Patched.txt` summarizes the TD burden of inherited LoFs in discovery cohort (SPARK+SSC). A related file `DiscovLoFTDBurdenExclKnown_Patched.txt` summarize the TD burden after excluding known ASD/NDD genes. Note: we used patched version throughout all analyses.
- `../../201912_RareVars/TDT_Patch/DiscovLoFTDBurdenChrX_Patched.txt` summarizes the TD burden of chrX LoFs. 
- `../../201912_RareVars/TDT_Patch/DiscovLoFCCBurden_Patched.txt` summarizes the pseudo CC burden.
- `../../201912_RareVars/TDT_Missense/DiscovDmisTDBurden.txt`summarize the TD burden of Dmis variants using different Dmis definitions.

The content of the above four sets of files were merged in `TDBurden/DiscovBurden.xlsx`. Similar sets of summary tables were generaeted for combined Discov+Rep cohorts, and the results were merged into `TDBurden/MetaBurden.xlsx` For discovery cohrot, we also evaluated burden of D-mis variants using different filtering criteria, and the results were merged to `TDBurden/DiscovBurdenDmis.xlsx`.

**Key observations:**

Burden were analyzed by using URLowAF LoFs with pExt>0.1. The results were visualized by comparing %transmission and over-transmission rate to affected vs unaffected offspring across different sets of genes. Similar to de novo variants, we found TD signals also concentrate in top two LOEUF deciles. But only 20% over-transmission to affected can be explained by known genes.

LoF TD was also found to be higher to low functionning ASD cases (defined as reported cognitive impairment in SPARK or Vineland<70 in SSC).

#### Gene set enrichment

We analyzed the enrichment of TD signals in the same sets of genes used for dnEnrich.

Results were summarized in `TDBurden/TDGeneSetEnrich.pptx` and `TDBurden/TDGeneSetEnrich.xlsx`

#### Gene table 

The final gene level summary for discovery cohort include the number LoF varinats (using different filters) in each gene for different sample groups and for variants in parents also includes transmission pattern. The initial results for autosome and chrX were separated into the following two files:

- `../../201912_RareVars/TDT_Patch/DiscovGeneTab_Patched.txt`

- `../../201912_RareVars/TDT_Patch/DiscovGeneTabChrX_Patched.txt `

For the purpose of prioritization, we only collected transmission information (by `TDBurden/scripts/prep_discov_genetab.pl`, and results are reformatted in TDBurden/TDTGenesDiscov.xlsx ). Pseudo case-controls will be merged with replication data in meta-analysis. 

We used the following criteria for prioritizing autosomal LoF-TD genes, for chrX only the top 1st gene was selected. Prioritized genes can given in `TDTGenesDiscov/Prioritized.txt`

```
LOEUFbin!="." & LOEUFbin<=1 & Arisk!="." & Arisk>=0.4
URLoFpExt0.1TDTStat>1 & AllLoFpExt0.1TDTStat>1
URLoFpExt0.1Par>0.5*URLoFpExt0.1Par
```

A total of 260 TD genes were prioritized (1 chrX), 15 overlap with genes ascertained from previous DenovoWEST analysis => 404 unique genes (13 chrX) carried forward for replication.

Note: privious analysis used slightly different criteria for prioritization but included all prioritized genes in this final version.

### Meta analysis

#### Gene list

We first create a full list of genes `Meta/CandGenes.txt` carried forward for meta-analysis (`Meta/scripts/create_meta_genelist.pl`). They include 404 genes in total (13 chrX) 138 of which are in known genes list. 159 of thoese genes are ascertained from de novo analysis, 260 from TDT analysis, 15 are shared in common. 379 genes are in LOEUF top 30%.

#### Burden

Overall burden was collected into `Meta/MetaBurden.xlsx`. We confirmed the main findings from discovery cohort, including association with ID.

#### CAFs of LoFs in controls

Evaluating the CAFs of different control sets on autosomal genes in top 30% LOEUF: `Meta/PopCAFs.pptx`. Genes CAFs that were fixed by manual curation were collected by `Meta/scripts/collect_gnomAD_fixedCAFs.pl`, and can be found in `Meta/PopCAFs.xlsx`.  Genes that have dramtic changes after curation were highlighted. Note: We do not show chromosome X, because for most genes LoF rate is too small. 

#### pExt thresholds

Before evaluating known gene burdens, we use de novo LoF enriched genes to calibrate gene specific pExt threshold. We have also removed variants that are curated nonLoFs by gnomAD.

We selected 96 known autosomal genes (plus MARK2)  with >=4 dnLoFs in ASD and NDD, in top 30% LOEUF,  then selected pExt threshold from 0.9, 0.5, and 0.1 that retain all LoFs with pExt>=0.1. Gene-specific thresholds can be found in `Meta/pExtThres.xlsx`, and CDFs of pExt in de novo and standing LoFs can be found in `Meta/pExtThres.pptx`. By applying those gene-specfic threshold and also removing curated nonLoFs, we found it further removes 18% standing LoFs that does not contribute to transmission disequilibrium.

#### Gene level meta-analysis

The initial results from meta- or mega- analysis were collected by the following script `Meta/scripts/collect_meta_results_raw.pl`, and can be found in `Meta/MetaGeneTab/RawResults.txt`. Raw results with fixed CAFs in population references are given in `Meta/MetaGeneTab/RawResults_Fixed.txt`. Note that we make CAF plots using LoF rates in those raw results. Initial raw and curated results are reformatted and can be found in `Meta/MetaGeneTab.xlsx`.

For dnLoF enriched genes, we applied variable pExt threshold and have curated gnomAD to remove Low conf LoF variants in population controls, results are collected into `Meta/MetaGeneTab/Curated.txt`. The curated and raw results are then merged, with priority given to curated results (`Meta/MetaGeneTab/MergeCuratedRaw.txt`). This will be the final results from meta and mega analysis, and reformatted to `Meta/MetaGeneTabFinal.xlsx`. We will draw final gene LoF distribution and case-control comparison plots using this data table.

For chrX, the results are collected into `Meta/MetaGeneTab/chrX.txt` and `Meta/MetaGeneTab/chrX_curated.txt` (the second file patched a few genes after removing FP or nonLoF calls in reference panel)

Note: In thoese tables, we included the number of de novo variants in full discov+rep cohort, or in all inhouse samples included for mega-analysis. For offspring samples included in TDT, we also collected number of inherited variants (from unaffected parents). Note that this differ slightly from the total number of inherited variants. Note also that the number of variants in de novo, inherited, and pseudo cases does not add up to the number in unreleated cases, because we may further remove related samples between TDT and pseudo case-controls (how many of them are related?)

```bash
# Coverages for autosomes are manually fixed into Excel.
table_intersect.pl -a Meta/MetaGeneTab.xlsx -a-sheet 2 -a-skip 1 -a-maxcol 26 -a-select HGNC -b ../../201912_RareVars/TDT_Patch/MetaGeneTab_Patched_stat.txt -b-select HGNC,SPARKOver15,gnomADexomeOver15,gnomADgenomeOver15,TopMedOver15 -b-fsep "\t" > Meta/MetaGeneTab/CovPatch.txt
table_intersect.pl -a Meta/MetaGeneTabFinal.xlsx -a-sheet 1 -a-skip 1 -a-maxcol 26 -a-select GeneID,HGNC -b ../../201912_RareVars/TDT_Patch/MetaGeneTab_Patched_stat.txt -b-select GeneID,HGNC,SPARKOver15,gnomADexomeOver15,gnomADgenomeOver15,TopMedOver15 -b-fsep "\t" > Meta/MetaGeneTab/CovPatchFinal.txt
# ChrX curated results were collected by following command:
perl Meta/scripts/collect_meta_results_chrX.pl ../../201912_RareVars/TDT_Patch/Curated/MetaGeneTabChrX_stat.txt Meta/MetaGeneTab/chrX_curated.txt
```

#### Gene level mega-analysis

The rationale for mega-analysis is to avoid the possible bias in meta-analysis that are driven by de novo enrichment especially when baseline mutation rate may deviate from presumed rate. We summarize the number of LoFs of three different types (de novo, inherited, and in pseudo controls), TDT signal, and rate comparison with three population references.  To make plots of mega-analysis, we reformat previous tables using `Meta/scripts/refmt_metatab4ccplot.pl`.

```bash
# First, collect raw results => Meta/MetaGeneTab/RawResults.txt
perl Meta/scripts/collect_meta_results_raw.pl
# Then, Collect results after fixing pop ref CAFs
# => Meta/MetaGeneTab/RawResults_Fixed.txt
perl Meta/scripts/collect_meta_results_raw.pl ../../201912_RareVars/TDT_Patch/Curated/MetaGeneTab_RawFixed_stat.txt Meta/MetaGeneTab/RawResults_Fixed.txt
# Also collect results after applying DNV guided thresholds
# => Meta/MetaGeneTab/Curated.txt
perl Meta/scripts/collect_meta_results_curated.pl
# For chrX, raw results
# => Meta/MetaGeneTab/chrX.txt
perl Meta/scripts/collect_meta_results_chrX.pl
# chrX: after fixing population CAFs
# => Meta/MetaGeneTab/chrX_curated.txt
perl Meta/scripts/collect_meta_results_chrX.pl ../../201912_RareVars/TDT_Patch/Curated/MetaGeneTabChrX_stat.txt Meta/MetaGeneTab/chrX_curated.txt
# Merge results by incoporating curated results into raw results
# We only keep DNV guided threshold genes with >=4 dnLoFs in ASD+DDD
# => Meta/MetaGeneTab/MergeCuratedRaw.txt
perl Meta/scripts/merge_raw_curated.pl
# Prepare tables for plots
perl Meta/scripts/refmt_metatab4ccplot.pl  Meta/MetaGeneTab/RawResults.txt
perl Meta/scripts/refmt_metatab4ccplot.pl  Meta/MetaGeneTab/MergeCuratedRaw.txt
perl Meta/scripts/refmt_metatab4ccplot.pl  Meta/MetaGeneTab/Curated.txt

```

To make plots, we also need to calculate sample sizes for three types of LoFs and genes that are curated.

```bash
# Note: for final pyramid plot, sample sizes need to be calculated.
perl Meta/scripts/calc_pyramid_sampsize.pl > Meta/MetaGeneTab/PlotSampSizes.txt
# Also find out list of genes that are curated
find ../../201912_RareVars/dat/PopRef/gnomADV211Exome/fix/  -name '*.gz' | while read vartab; do basename $vartab | sed -s 's/.gz//'; done | awk 'BEGIN{print "HGNC"}{print $0}' > Meta/MetaGeneTab/ExomeFixed.txt
find ../../201912_RareVars/dat/PopRef/gnomADgenomeV31WGSCoding/fix/  -name '*.gz' | while read vartab; do basename $vartab | sed -s 's/.gz//'; done | awk 'BEGIN{print "HGNC"}{print $0}' > Meta/MetaGeneTab/GenomeFixed.txt
# Making plots
Rscript Meta/MetaGeneTab.R
Rscript Meta/MetaGeneTab_New.R
Rscript Meta/MetaGeneTabLong.R
Rscript Meta/MetaGeneTab_chrX.R
```

#### Data

```bash
# First output TDT samples
perl Data/scripts/refmt_TDTsamptabs.pl
# Then reformat TDT vars
perl Data/scripts/refmt_TDTvartabs.pl
# Output Pseudo CC samples
perl Data/scripts/refmt_PseudoCCsamptabs.pl
# Reformat Pseudo CC vars
perl Data/scripts/refmt_PseudoCCvartabs.pl
# Mega cases list
perl Data/scripts/refmt_Megasamtabs.pl
# Mega LoFs
perl Data/scripts/refmt_MegaCasevartabs.pl

```

