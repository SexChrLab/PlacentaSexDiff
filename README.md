# PlacentaSexDiff
Sex differences in human placentas

#### Allele specific expression analysis
* Data: 
  - 11 female placentas
  - GATK Haplotype single sample genotype calling
  - Processed RNAseq bam file. For each sample, there are 2 placenta samples for 2 different location

* Steps:
1. Run ASEReadcounter:
  - See the snakemake file `exome.snakefile`.
2. Filter the results from ASEReadCounter to contain just the heterozygous sites
3. Plot histogram of Alternate allele count / Total count ratio for autosomes and X chromosome
4. Compute directionality:
  - If Alternate allele count / Total count ratio < 0.5 for both locations, annotate "Same"
  - If Alternate allele count / Total count ratio < 0.5 for one location but > 0.5 for the other location, and vice versa, annotate "Diff".
5. Plot scatter plot of Alteranate allele count / Total count ratio for 2 locations.
