# PlacentaSexDiff
Sex differences in human placentas


## Sex differences in placenta expression 
#### Step 1: Process RNA sequence data
- FastQC, MultiQC, and bamstats for visualizing quality. bbduk to trim reads and remove adaptors. 
- Reads aligned to a GRCh38 sex chromosomome complement reference genome. See https://github.com/SexChrLab/XY_RNAseq for more details on aligning to a sex chromosomome complement reference genome
- Use the snakemake file `rna.snakefile` that includes the following steps:
 - trimming, alignmnet, obtaining counts, and generating bam stats. 

#### Step 2: Run Limma/Voom 
- Use the r file `R_limmaVooom_placentaSexDiff.r` that includes the following steps:
  - You will need a counts table generated from step 1 FeatureCounts and a pheno type file
  - run limma/voom


## Heterogeneity
#### Step 1: Process exome sequence data
- Previously, we processed batch 1 and batch 2 separately. However, we thought that genotype calling will have more power if we combine all of the exome sequence data from batch 1 and batch 2 together.
- Previously, we already mapped to sex-specific reference genome and produced the gVCF file. Therefore, here we need to combine all of the gVCF files from both batches and run GenotypeGVCF.
- Use the snakemake file `exome.snakefile` that includes the following steps:
 - joint-genotyping across all 30 samples

#### Step 2: Run GATK ASEReadCounter
- Use the snakemake file `asereadcounter.snakefile` that includes the following steps:
  - make symlink for RNAseq bam files (processed by HISAT)
  - make symlink for VCF
  - run GATK ASEReadCounter
- Note that when I run ASEReadCounter, I am using the VCF file immmediately after GATK's joint genotyping (i.e. before any filtering).
- Run GATK ASEReadCounter for 3 threshold of minDepth: 10, 20, and 50

#### Step 3: Process outputs from ASEReadCounter
- Because in Step 2, I use the VCF file immmediately after GATK's joint genotyping, in this step, I filter the outputs from ASEReadCounter to contains only biallelic and heterozygous sites
  - Use the Python script `obtain_vqsrfiltered_biallelic_hets_sites.py`.
- For some positions, ref allele and alternate alleles are the same and are not filtered by previous filtering. So I write a Python script to remove these sites.
  - Use the Python script `remove_same_refAlt.py`.
- Use the snakemake file `post_asereadcounter_processing.snakefile`.
