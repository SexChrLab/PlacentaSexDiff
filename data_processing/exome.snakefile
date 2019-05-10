import os

configfile: "exome_config.json"

# Tool paths:
fastqc_path = "fastqc"
multiqc_path = "multiqc"
bwa_path = "bwa"
samtools_path = "samtools"
bbduksh_path = "bbduk.sh"
picard_path = "picard"
gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk"
gatk3_path = "/home/tphung3/softwares/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar"

#
REF_TYPE = ["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]

# Directory
fastq_directory = "fastq_files/"

rule all:
    # input:
    #     expand("asereadcounter/{female_sample}/{female_sample}_placenta_2.csv", female_sample = config["females_placenta"])
    # input:
    #     expand("asereadcounter/{female_sample}/{female_sample}_placenta_1.csv", female_sample = config["females_placenta"])
    # input:
    #     expand("rnaseq_bam/{sample_name}/{sample_name}_placenta_1.bam.bai", sample_name = config["females_placenta"]),
    #     expand("rnaseq_bam/{sample_name}/{sample_name}_placenta_2.bam.bai", sample_name = config["females_placenta"]),
    input:
        expand("genotyped_vcfs/{female_sample}/{female_sample}.gatk.called.raw.vcf.gz", female_sample = config["females"])
    input:
        expand("gvcfs/{chrm}/{chrm}.{female_sample}.merged.g.vcf.gz", chrm=config["chr_females"], female_sample =config["females"])
    input:
        expand("processed_bams/{female_sample}.GRCh38_Ymasked.sorted.mkdup.bam.bai", female_sample = config["females"])
    input:
        expand("processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.mkdup.bam.bai", male_sample = config["males"])
    input:
        expand("processed_bams/{female_sample}.GRCh38_Ymasked.sorted.mkdup.bam", female_sample = config["females"]),
        expand("stats/{female_sample}.GRCh38_Ymasked.picard_mkdup_metrics.txt", female_sample = config["females"])
    input:
        expand("processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.mkdup.bam", male_sample = config["males"]),
        expand("stats/{male_sample}.GRCh38_minusYPARs.picard_mkdup_metrics.txt", male_sample = config["males"])
    input:
        expand("processed_bams/{female_sample}.GRCh38_Ymasked.sorted.bam.bai", female_sample = config["females"])
    input:
        expand("processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.bam.bai", male_sample = config["males"])
    input:
        expand("processed_bams/{female_sample}.GRCh38_Ymasked.sorted.bam", female_sample = config["females"])
    input:
        expand("processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.bam", male_sample = config["males"])
    # input:
    # 	"multiqc_results/multiqc_report.html"
    input:
        expand("trimmed_fastqs/{sample_name}_trimmed_R1.fastq.gz", sample_name = config["sample_names"]),
        expand("trimmed_fastqs/{sample_name}_trimmed_R2.fastq.gz", sample_name = config["sample_names"])
    input:
        expand("fastq_files/{sample_name}_R1.fastq.gz", sample_name = config["sample_names"])
    input:
        expand("refs/{ref_type}.fa.fai", ref_type = REF_TYPE),
        expand("refs/{ref_type}.fa.amb", ref_type = REF_TYPE),
        expand("refs/{ref_type}.dict", ref_type = REF_TYPE)



rule prep_refs_mk_sy_ln:
    input:
        ref = lambda wildcards: config["genome_paths"][wildcards.ref_type]
    output:
        ref_sy_ln = "refs/{ref_type}.fa"
    shell:
        """
        ln -s {input.ref} {output.ref_sy_ln}
        """

rule prep_refs:
    input:
        ref = "refs/{ref_type}.fa"
    output:
        fai = "refs/{ref_type}.fa.fai",
        amb = "refs/{ref_type}.fa.amb",
        dict = "refs/{ref_type}.dict"
    params:
        samtools = samtools_path,
        bwa = bwa_path
    run:
        # faidx
        shell("{params.samtools} faidx {input.ref}")

        # .dict
        shell("{params.samtools} dict -o {output.dict} {input.ref}")

        # bwa
        shell("{params.bwa} index {input.ref}")


rule mk_sy_ln_fastqs:
    input:
        original_R1 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq1"],
        original_R2 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq2"]
    output:
        R1_out = "fastq_files/{sample_name}_R1.fastq.gz",
        R2_out = "fastq_files/{sample_name}_R2.fastq.gz"
    shell:
        """
        ln -s {input.original_R1} {output.R1_out};
        ln -s {input.original_R2} {output.R2_out}
        """

rule trim_adapters_paired_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample_name]["fq1_sy"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample_name]["fq2_sy"])
	output:
		out_fq1 = "trimmed_fastqs/{sample_name}_trimmed_R1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample_name}_trimmed_R2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=/mnt/storage/SAYRES/REFERENCE_GENOMES/adapters/adapter_sequence.fa "
		"qtrim=rl trimq=30 minlen=75 maq=20"

##############################################################################
# The fastqc and multiqc rules are not working for now. Will need to come back
##############################################################################

# rule fastqc_analysis:
#     input:
#         fq1 = lambda wildcards: os.path.join(
# 			fastq_directory, config[wildcards.sample_name]["fq1_sy"]),
#         fq2 = lambda wildcards: os.path.join(
# 			fastq_directory, config[wildcards.sample_name]["fq2_sy"])
#     output:
#         fq1_fastqc = "fastqc_results/{sample_name}.R1_fastqc.html",
#         fq2_fastqc = "fastqc_results/{sample_name}.R2_fastqc.html"
#
#     params:
#         fastqc = fastqc_path
#
#     shell:
#         """
#         PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input.fq1};
#         PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input.fq2}
#         """

# rule fastqc_analysis:
# 	input:
# 		"fastq_files/{sample_name}_{read}.fastq.gz"
# 	output:
# 		"fastqc_results/{sample_name}.{read}_fastqc.html"
# 	params:
# 		fastqc = fastqc_path
# 	shell:
# 		"PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input}"
#
# rule multiqc_analysis:
# 	input:
# 		expand(
# 			"fastqc_results/{sample_name}.{read}_fastqc.html",
# 			sample_name=config["sample_names"],
# 			read=["R1", "R2"])
# 	output:
# 		"multiqc_results/multiqc_report.html"
# 	params:
# 		multiqc = multiqc_path
# 	shell:
# 		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
# 		"{params.multiqc} --interactive -f "
# 		"-o multiqc_results fastqc_results"

rule map_to_sex_specific_refs_males:
    input:
        fq1 = "trimmed_fastqs/{male_sample}_trimmed_R1.fastq.gz",
        fq2 = "trimmed_fastqs/{male_sample}_trimmed_R2.fastq.gz",
        ref = "refs/Ref_GRCh38_Y_PARsMasked.fa",
        fai = "refs/Ref_GRCh38_Y_HardMasked.fa.fai"
    output:
        "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.bam"
    params:
        id = lambda wildcards: config[wildcards.male_sample]["ID"],
        sm = lambda wildcards: config[wildcards.male_sample]["SM"],
        lb = lambda wildcards: config[wildcards.male_sample]["LB"],
        pu = lambda wildcards: config[wildcards.male_sample]["PU"],
        pl = lambda wildcards: config[wildcards.male_sample]["PL"],
        bwa = bwa_path,
        samtools = samtools_path,
        threads = 4
    threads: 4
    priority: 100
    shell:
        " {params.bwa} mem -R "
        "'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
        "{input.ref} {input.fq1} {input.fq2}"
        "| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
        "-O bam -o {output}"

rule map_to_sex_specific_refs_females:
    input:
        fq1 = "trimmed_fastqs/{female_sample}_trimmed_R1.fastq.gz",
        fq2 = "trimmed_fastqs/{female_sample}_trimmed_R2.fastq.gz",
        ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
        fai = "refs/Ref_GRCh38_Y_HardMasked.fa.fai"
    output:
        "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.bam"
    params:
        id = lambda wildcards: config[wildcards.female_sample]["ID"],
        sm = lambda wildcards: config[wildcards.female_sample]["SM"],
        lb = lambda wildcards: config[wildcards.female_sample]["LB"],
        pu = lambda wildcards: config[wildcards.female_sample]["PU"],
        pl = lambda wildcards: config[wildcards.female_sample]["PL"],
        bwa = bwa_path,
        samtools = samtools_path,
        threads = 4
    threads: 4
    priority: 100
    shell:
        " {params.bwa} mem -R "
        "'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
        "{input.ref} {input.fq1} {input.fq2}"
        "| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
        "-O bam -o {output}"

rule index_bam_males:
    input:
        "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.bam"
    output:
        "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule index_bam_females:
    input:
        "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.bam"
    output:
        "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule picard_mkdups_males:
    input:
        bam = "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.bam",
        bai = "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.bam.bai"
    output:
        bam = "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.mkdup.bam",
        metrics = "stats/{male_sample}.GRCh38_minusYPARs.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule picard_mkdups_females:
    input:
        bam = "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.bam",
        bai = "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.bam.bai"
    output:
        bam = "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.mkdup.bam",
        metrics = "stats/{female_sample}.GRCh38_Ymasked.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule index_mkdup_bam_males:
    input:
        "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.mkdup.bam"
    output:
        "processed_bams/{male_sample}.GRCh38_minusYPARs.sorted.mkdup.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule index_mkdup_bam_females:
    input:
        "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.mkdup.bam"
    output:
        "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.mkdup.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

# GATK for 12 female sample
rule gatk_gvcf_females:
	input:
		ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
		bam = "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.mkdup.bam",
		bai = "processed_bams/{female_sample}.GRCh38_Ymasked.sorted.mkdup.bam.bai"
	output:
		"gvcfs/{female_sample}/{female_sample}.g.vcf.gz"
	params:
		gatk = gatk_path,
	shell:
		"{params.gatk} "
		"HaplotypeCaller -R {input.ref} -I {input.bam} "
		"--emit-ref-confidence GVCF --output {output}"

rule gatk_genotypegvcf_females:
    input:
        ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
        gvcf = "gvcfs/{female_sample}/{female_sample}.g.vcf.gz"
    output:
        "genotyped_vcfs/{female_sample}/{female_sample}.gatk.called.raw.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """{params.gatk} --java-options "-Xmx10g" """
        """GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} """

###############################
# Filter with VQSR
# This currently is not working
###############################

# rule gatk_variantrecalibrator:
#     input:
#         ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
#         vcf = "genotyped_vcfs/{female_sample}/{female_sample}.gatk.called.raw.vcf.gz",
#         hapmap = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/hapmap_3.3.hg38.vcf.gz",
#         omni = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/1000G_omni2.5.hg38.vcf.gz",
#         thousandG = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
#         dbsnp = "/scratch/tphung3/PopulationReferenceAlignment/analyses/compare_hardfilter_vqsr/scripts/dbsnp_138.hg38.vcf.gz"
#     output:
#         recal = "vqsr/{female_sample}/{female_sample}_output.recal",
#         tranches = "vqsr/{female_sample}/{female_sample}_output.tranches"
#     params:
#         gatk = gatk_path
#     shell:
#         """{params.gatk} --java-options "-Xmx16g" VariantRecalibrator """
#         """-R {input.ref} -V {input.vcf}  """
#         """--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} """
#         """--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} """
#         """--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} """
#         """--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} """
#         """-an QD -an ReadPosRankSum -an FS -an SOR """
#         """-mode SNP """
#         """-O {output.recal} """
#         """--tranches-file {output.tranches} """
#
# rule gatk_applyvqsr:
#     input:
#         ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
#         vcf = "genotyped_vcfs/{female_sample}/{female_sample}.gatk.called.raw.vcf.gz",
#         tranches = "vqsr/{female_sample}/{female_sample}_output.tranches",
#         recal = "vqsr/{female_sample}/{female_sample}_output.recal"
#     output:
#         "vqsr/{female_sample}/{female_sample}.gatk.called.vqsr.vcf.gz"
#     params:
#         gatk = gatk_path
#     shell:
#         """{params.gatk} --java-options "-Xmx16g" ApplyVQSR """
#         """-R {input.ref} """
#         """-V {input.vcf} """
#         """-O {output} """
#         """--truth-sensitivity-filter-level 99.0 """
#         """--tranches-file {input.tranches} """
#         """--recal-file {input.recal} """
#         """-mode SNP """
#
# rule gatk_selectvariants:
#     input:
#         ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
#         vcf = "vqsr/{female_sample}/{female_sample}.gatk.called.vqsr.vcf.gz"
#     output:
#         "vqsr/{female_sample}/{female_sample}.gatk.called.vqsr.sv.vcf.gz"
#     params:
#         gatk = gatk_path
#     shell:
#         """{params.gatk} --java-options "-Xmx16g" SelectVariants """
#         """-R {input.ref} """
#         """-V {input.vcf} """
#         """--exclude-filtered """
#         """-O {output} """

#####################################################################################
# Run ASEReadCounter.
# However, the bam files used here were from a different reference. This is archived.
# Use a different snakemake file for the ASEReadCounter step.
#####################################################################################

# # Convert the RNA seq bam files to have chr1 format instead of 1
# rule convert_bam_chr_fmt:
#     input:
#         placenta_1_bam = lambda wildcards: config[wildcards.sample_name]["placenta_1_bam_path"],
#         placenta_2_bam = lambda wildcards: config[wildcards.sample_name]["placenta_2_bam_path"]
#     output:
#         placenta_1_bam_fmt = "rnaseq_bam/{sample_name}/{sample_name}_placenta_1.bam",
#         placenta_2_bam_fmt = "rnaseq_bam/{sample_name}/{sample_name}_placenta_2.bam"
#     params:
#         samtools = samtools_path
#     shell:
#         """
#         {params.samtools} view -H {input.placenta_1_bam} | sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:MT/SN:chrM/' | {params.samtools} reheader - {input.placenta_1_bam} > {output.placenta_1_bam_fmt};
#         {params.samtools} view -H {input.placenta_2_bam} | sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:MT/SN:chrM/' | {params.samtools} reheader - {input.placenta_2_bam} > {output.placenta_2_bam_fmt}
#         """
#
# rule index_bam_chr_fmt:
#     input:
#         placenta_1_bam_fmt = "rnaseq_bam/{sample_name}/{sample_name}_placenta_1.bam",
#         placenta_2_bam_fmt = "rnaseq_bam/{sample_name}/{sample_name}_placenta_2.bam"
#     output:
#         placenta_1_bam_fmt_bai = "rnaseq_bam/{sample_name}/{sample_name}_placenta_1.bam.bai",
#         placenta_2_bam_fmt_bai = "rnaseq_bam/{sample_name}/{sample_name}_placenta_2.bam.bai"
#     params:
#         samtools = samtools_path
#     shell:
#         """
#         {params.samtools} index {input.placenta_1_bam_fmt};
#         {params.samtools} index {input.placenta_2_bam_fmt}
#         """
#
# # ASEReadCounter on the no filter VCF file
# rule gatk_asereadcounter_placenta1:
#     input:
#         ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
#         bam = "rnaseq_bam/{female_sample}/{female_sample}_placenta_1.bam",
#         sites = "genotyped_vcfs/{female_sample}/{female_sample}.gatk.called.raw.vcf.gz"
#     output:
#         "asereadcounter/{female_sample}/{female_sample}_placenta_1.csv"
#     params:
#         gatk = gatk3_path
#     shell:
#         """java -jar {params.gatk} """
#         """-R {input.ref} """
#         """-T ASEReadCounter """
#         """-o {output} """
#         """-I {input.bam} """
#         """-sites {input.sites} """
#         """-U ALLOW_N_CIGAR_READS """
#         """-minDepth 10 """
#         """--minMappingQuality 10 """
#         """--minBaseQuality 2 """
#         """-drf DuplicateRead """
#
# rule gatk_asereadcounter_placenta2:
#     input:
#         ref = "refs/Ref_GRCh38_Y_HardMasked.fa",
#         bam = "rnaseq_bam/{female_sample}/{female_sample}_placenta_2.bam",
#         sites = "genotyped_vcfs/{female_sample}/{female_sample}.gatk.called.raw.vcf.gz"
#     output:
#         "asereadcounter/{female_sample}/{female_sample}_placenta_2.csv"
#     params:
#         gatk = gatk3_path
#     shell:
#         """java -jar {params.gatk} """
#         """-R {input.ref} """
#         """-T ASEReadCounter """
#         """-o {output} """
#         """-I {input.bam} """
#         """-sites {input.sites} """
#         """-U ALLOW_N_CIGAR_READS """
#         """-minDepth 10 """
#         """--minMappingQuality 10 """
#         """--minBaseQuality 2 """
#         """-drf DuplicateRead """
