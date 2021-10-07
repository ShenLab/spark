### Known genes

Known ASD/NDD genes are used to look for likely pathogenic variants and to evaluate "residual" burden of de novo and transmitted variants. We may consider allelic requirement and mutation consequences for different types of analyses. For example, when analyzing residual de novo/TDT burden, we considered all genes that are likely dominant. But in pratice allelic requirement is not used, because in many genes, full allelic spectrum of causative mutations have not been established. 

#### Known dominant NDD genes

We collected dominant NDD genes from DDG2P database (2020-02). Two list is prepared, the stringent list only include genes that are labeled as *confirmed* or *probable*.

```bash
# Known NDD genes (591 uniq coding): no filter on allelic requirement!
flatdb.pl --in ~/Dropbox/Data/AnnoDB/DDG2P/DDG2P_26_2_2020_wGeneID.tsv -f 'DDDAllelic=~"mono|x-link|hemi" &  (DDDOrgan=~"Brain" | DDDOrgan=="")' -fsep "\t" --alias gencodev19_id:GeneID,gene_symbol:HGNC,allelic_requirement:DDDAllelic,mutation_consequence:DDDMutCons,organ_specificity_list:DDDOrgan,DDD_category:DDDCategory --select GeneID,HGNC,DDDCategory,DDDAllelic,DDDMutCons,DDDOrgan | table_intersect.pl  -a - -a-fsep "\t" -b ~/Dropbox/Data/GeneMetrics/SFARI_release2019-11-13.txt -b-fsep "\t" -b-select GencodeV19ID:GeneID,Category:SFARICategory,GeneScore:SFARIScore -with-loj | table_intersect.pl -a - -a-fsep "\t" -b ~/Dropbox/Data/GeneMetrics/SFARI_release2019-08-29.txt -b-fsep "\t" -b-select GencodeV19ID:GeneID,GeneScore:SFARIScore2 -with-loj  > Known/DomNDDGenes.txt
# A more stringent list (508 uniq coding), including only "confirmed/probable" 
flatdb.pl --in Known/DomNDDGenes.txt -f 'DDDCategory=~"confirmed|probable"' -fsep "\t" --out Known/DomNDDGenesStringent.txt
```

#### ASD genome-wide signficant genes

They include a total of 100 genes that are enriched with de novo variants by p<1e-3 in autism discovery cohort (17K trios) and reached genome-wide signfiicant in full autism cohort (23K tiros) or in autism + NDD trios (54K trios). The genes list are given in `AutismGWSGenes.txt` (by **Known/scripts/collect_autism_gwsgenes.pl**). This list include novel autism/NDD genes like *MARK2*, *PRKAR1B*. We manually removed two genes *SMAD6* and *TNPO3* from this list because the formal has high LoF rate in population, and evidence for the latter from de novo+TDT+CC combined meta-analysis is ambiguous .

From this gene set, we noted that:

* Most ASD genes are also nominally significant in NDD data alone, except *SMARCC2* and *BRSK2*.
* All genes that reached genome-wide signficance (GWS) in autism discovery cohort remain GWS in autism or autism+NDD meta-analysis.

* The contribution of LoF and Dmis to the disease are generally concordant between ASD and NDD. A few exceptions are observed in channel genes which are consistent with known genotype-phenotype association (More on this point later in analyzing top genes).

#### SPARK gene list

We compared the above two sets of genes with SPARK genes (2020Jul, `Known/SPARKGenes.txt`, by Known/scripts/compare_SPARK_genes.pl). Out of 157 SPARK genes, 24 were not included in the above two list, including 10 known NDD genes with biallelic inheritance, and 14 other genes curated from autism literature (Curated in `Known/SPARKGenes.xlsx`).

#### Final gene list

We will define two sets of "known" autism/NDD genes for de novo/TDT analysis:

* The **stringent** set (`Known/GeneList_Strict.txt`, n=511) includes: stringent dominant NDD genes + autism genes reached genome-wide signficance (excluding *TNPO3* and *SMAD6* that have inconclusive evidence) 

* The **comprehensive** set (`Known/GeneList_Full.txt`, n=618) include: all dominant NDD genes + all autism genes reached genome-wide signficance  + all dominant SPARK genes + we also not include additional genes ranked by TADA that have SFARI score 1 or 2 (Note: new genes identified this study were also not included as known.)

Also generate a list of  61 autism GWS genes (`Known/GeneList_ASDGWS.txt` filtered by pASDMeta<1.3e-6). 

### Functional sets

We collected following functional gene sets from the literature for gene set analysis.

- **Transcriptome**: Brain-enriched, Neuron-specific, Inhibitory neuron-specific, Excitatory neuron-specific, Co-expression modules M2/M3
- **Regulome**: CELF4, CHD8, FMRP, RBFOX1/3, RBFOX2 target genes
- **Proteome**: SynaptomeDB

### Genetic sets

We prepared two classes of gene sets based on genetic metrics.

- Genetic evidence: including genes enriched with de novo variants in ASD or NDD (using nomial p-value cutoffs), and genes enriched 

- Constrained genes: based on ExAC sHet or gnomAD LOEUF

For gene set enrichment analysis. Both functional and genetic sets are further restricted to a background set of **5804** genes (5754 are callable and not excluded by blacklist) constrained genes defined as pLI>=0.5 or top 20% LOEUF. 

We reformatted gene sets into membership matrix format and count the number of genes in used in dnEnrich and TDenrich analysis (by `FuncGenet/scripts/collect_genesets.pl`). Results can be found in `FuncGenet/GeneSets.xlsx` spreadsheet. The pairwise overlap between gene sets were also visualized (by `FuncGenet/GeneSets_PairwiseOverlaps.R`) and figure can be found in `FuncGenet/GeneSets.pptx`. 

In burden plots, we also used ExAC pLI or gnomAD LOEUF deciles to stratify genes. The number of genes in each set will be shown in the burden plots. So we also tally the number of genes (by `FuncGenet/scripts/tally_burden_genesets.pl`), results can be found in `FuncGenet/GeneSets4Burden_NumGenesInSet.txt`.

