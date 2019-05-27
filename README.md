# PlacentaSexDiff
Sex differences in human placentas

### Explore heterogeneity in allele specific expression

**Pilot analysis on placenta 1 and placenta 2 from sample OBG0044**

* Exome:
  - GATK Haplotype single sample genotype calling
* RNAseq:
  - Processed by HISAT
  
* Both Exome and RNAseq were mapped to the same refrence 

* Steps:
1. Run GATK ASEReadCounter:
* Use the snakemake file `asereadcounter.snakefile`.

##### Per variant
* Steps:
1. Run ASEReadcounter:
  - See the snakemake file `exome.snakefile`.
2. Filter the results from ASEReadCounter to contain just the heterozygous sites
  - Use the python script `filter_for_hets.py`
3. Plot histogram of Alternate allele count / Total count ratio for autosomes and X chromosome & merge placenta 1 and placenta 2 for directionality analysis.
4. Compute directionality:
  - If Alternate allele count / Total count ratio < 0.5 for both locations, annotate "Same"
  - If Alternate allele count / Total count ratio < 0.5 for one location but > 0.5 for the other location, and vice versa, annotate "Diff".
  
  ```
  python compute_bias_directionality.py OBG0044_placenta_1_2_hets.csv OBG0044_placenta_1_2_hets_directionality.csv
  ```
5. Plot scatter plot of Alteranate allele count / Total count ratio for 2 locations.

##### Per gene
* Steps:
1. Use VEP to annotate the gene for each variant (need to set this up on the cluster).
2. Use `find_num_snps_per_gene.py` to obtain the alternate allele count and total allele count for each gene on the autosomes and chrX. Right now, it's summing up the alternate allele count and total count from all of the variants for that gene. 
3. Use the Rscript `analyze_ase.R` to merge placenta 1 and placenta 2. 
4. Use the script `comppute_bias_directionality_gene.py` to find the direction in bias. 
5. Plot using the Rscript `analyze_ase.R`

#### PHASING
##### Phase on the whole chrX
* Use the script `phasing.py`. 

```
python phasing.py OBG0044_placenta_1_hets_totalcountgreater20.csv OBG0044_placenta_2_hets_totalcountgreater20.csv placenta_1_hap_1_ratio.txt placenta_2_hap_1_ratio.txt

python phasing.py OBG0044_placenta_1_hets_totalcountgreater20.csv OBG0044_placenta_2_hets_totalcountgreater20.csv placenta_1_hap_2_ratio.txt placenta_2_hap_2_ratio.txt
```

##### Phase on genes where each gene has at least 2 variants
* Use the script `phasing_genes_more_than_2_variants.py`.
* Note: need to change line 57 and 62 for different haplotypes

```
python phasing_genes_more_than_2_variants.py OBG0044_placenta_1_hets_totalcountgreater20.csv OBG0044_placenta_2_hets_totalcountgr
eater20.csv OBG0044.gatk.called.raw_vep.vcf TEST_GENE_placenta_1_hap_1.txt TEST_GENE_placenta_2_hap_1.txt"

python phasing_genes_more_than_2_variants.py OBG0044_placenta_1_hets_totalcountgreater20.csv OBG0044_placenta_2_hets_totalcountgr
eater20.csv OBG0044.gatk.called.raw_vep.vcf TEST_GENE_placenta_1_hap_2.txt TEST_GENE_placenta_2_hap_2.txt"
```

## PART A. For variants and genes that show allele specific expression, is the haplotype that are expressed the same in both extraction sites from the same placenta? 

* Data: 12 placenta samples from 12 individuals. Exome on 12 placentas. RNA-seq on 2 extraction sites on each placenta. 

### Whole exome processing
* Joint call using GATK4 across 12 individuals. 
  - Subset the VCF file into 2 files: chrA and chrX
  
  ```
  bcftools view batch_1.gatk.called.raw.females.vcf.gz --regions chrX > batch_1.gatk.called.raw.females.chrX.vcf
  bcftools view batch_1.gatk.called.raw.females.vcf.gz --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 > batch_1.gatk.called.raw.females.chrA.vcf
  ```
  
* Filter using VQSR
* After filteirng using VQSR, filter for biallelic variants only. Then, subset for each individual. After this step, there would be variants that are homozygous reference or homozygous alternate here. Therefore, run a GATK command to obtain the heterozygous sites. 

* Examine annotations such as DP and MQ on VCF file 
```
python ~/softwares/tanya_repos/vcfhelper/extract_stats_from_vcf.py QD FS SOR MQ MQRankSum ReadPosRankSum --vcf batch_1.gatk.called.raw.females.chrA.vcf --outfile ../analysis/post_gatk_call_processing_exome/annotations/chrA_prefiltering_annotations.txt
```

### Run GATK ASEReadCounter
* RNA-seq bam files are from HISAT
* Variant files are after the GATK step of joint calling.
**NOTES: I tried running GATK ASEReadCounter using the variant files after filtering for biallelic sites, subset per individual, and filtering for hets sites. However, the command kept failing, saying that the contigs for the bam files and the refenrece do not match. However, using the variants files before these processing steps work.**

### Filter results from GATK ASEReadCounter for hets
