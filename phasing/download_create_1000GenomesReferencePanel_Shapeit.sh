#!/bin/bash

date

cd /mnt/storage/SAYRES/placenta_YPOPS/07_ASE/phase_with_1000Genomes/1000Genomes_reference_panel

# required packages
module load bamtools/2.4.0 
module load java/latest
module load bcftools/1.4.0 
module load vcftools/0.1.12b
module load shapeit/v2.r837

#--------
#   1. download 1000Genomes VCF lifted over to GRCh38 
#       Located here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
#--------

# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi


wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr3.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr3.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr4.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr4.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr7.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr7.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr9.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr9.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr15.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr15.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr16.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr16.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr18.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr18.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr19.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr19.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

#--------
#   2. Convert 1000Genomes.GRCh38.vcf to .hap and .sample 
#
#--------

# --hapsample prefix or hap-file,sample-file
# convert from VCF to hap/sample format used by IMPUTE2 and SHAPEIT. The columns of .hap file begin with ID,RSID,POS,REF,ALT. In order to prevent strand swaps, the program uses IDs of the form "CHROM:POS_REF_ALT".

bcftools convert --haplegendsample chr1 ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr2 ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr3 ALL.chr3.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr4 ALL.chr4.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr5 ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr6 ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr7 ALL.chr7.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr8 ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr9 ALL.chr9.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr10 ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr11 ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr12 ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr13 ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr14 ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr15 ALL.chr15.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr16 ALL.chr16.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr17 ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr18 ALL.chr18.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr19 ALL.chr19.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr20 ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr21 ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chr22 ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
bcftools convert --haplegendsample chrX ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

#--------
date
