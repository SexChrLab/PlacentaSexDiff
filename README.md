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
