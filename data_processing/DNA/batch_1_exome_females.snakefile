configfile: "batch_1_exome_config.json"

# conda environment: placenta

# Tool paths:
fastqc_path = "fastqc"
multiqc_path = "multiqc"
bwa_path = "bwa"
samtools_path = "samtools"
bbduksh_path = "bbduk.sh"
picard_path = "picard"
gatk_path = "/home/tphung3/softwares/gatk-4.1.0.0/gatk" #change the path here for GATK4
gatk3_path = "/home/tphung3/softwares/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar" #change the path here for GATK3
bcftools_path = "bcftools"

# Directory
fastq_directory = "fastq_files/"

rule all:
    input:
        expand("processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.mkdup.bam.bai", female_sample=config["females"])

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

# ----------------------------------
# FASTQC and MULTIQC before trimming
# ----------------------------------

rule fastqc_analysis:
    input:
        fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample_name]["fq1_sy"]),
        fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample_name]["fq2_sy"])
    output:
        fq1_fastqc = "fastqc_results/{sample_name}.R1_fastqc.html",
        fq2_fastqc = "fastqc_results/{sample_name}.R2_fastqc.html"

    params:
        fastqc = fastqc_path

    shell:
        """
        PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input.fq1};
        PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input.fq2}
        """

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample_name}.{read}_fastqc.html",
			sample_name=config["sample_names"],
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

# ----------------------------------
# FASTQC and MULTIQC after trimming
# ----------------------------------
rule fastqc_analysis_post_trimming:
    input:
        fq1 = "/data/storage/SAYRES/placenta_YPOPS/08_exome/trimmed_fastqs/{sample_name}_trimmed_R1.fastq.gz",
        fq2 = "/data/storage/SAYRES/placenta_YPOPS/08_exome/trimmed_fastqs/{sample_name}_trimmed_R2.fastq.gz"
    output:
        fq1_fastqc = "fastqc_results/{sample_name}_trimmed_R1_fastqc.html",
        fq2_fastqc = "fastqc_results/{sample_name}_trimmed_R2_fastqc.html"

    params:
        fastqc = fastqc_path

    shell:
        """
        PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input.fq1};
        PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input.fq2}
        """

# move the trimmed fastqc results files to the directory fastqc_results_trimmed
rule multiqc_analysis_trimming:
	input:
		expand(
			"fastqc_results/{sample_name}_trimmed_{read}_fastqc.zip",
			sample_name=config["sample_names"],
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

# ----------------------------------
# Exome processing for females samples
# ----------------------------------
rule prep_refs_mk_sy_ln:
    input:
        fa = config["XX_ref_fa"],
        fai = config["XX_ref_fai"]
    output:
        fa = "refs/GRCh38.p12.genome.XXonly.fa",
        fai = "refs/GRCh38.p12.genome.XXonly.fa.fai"
    shell:
        """
        ln -s {input.fa} {output.fa};
        ln -s {input.fai} {output.fai}
        """

rule prep_refs:
    input:
        ref = "refs/GRCh38.p12.genome.XXonly.fa"
    output:
        amb = "refs/GRCh38.p12.genome.XXonly.fa.amb"
    params:
        bwa = bwa_path
    run:
        # bwa
        shell("{params.bwa} index {input.ref}")

rule map_to_sex_specific_refs_females:
    input:
        fq1 = "trimmed_fastqs/{female_sample}_trimmed_R1.fastq.gz",
        fq2 = "trimmed_fastqs/{female_sample}_trimmed_R2.fastq.gz",
        ref = "refs/GRCh38.p12.genome.XXonly.fa",
        fai = "refs/GRCh38.p12.genome.XXonly.fa.fai"
    output:
        "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.bam"
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

rule index_bam_females:
    input:
        "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.bam"
    output:
        "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule picard_mkdups_females:
    input:
        bam = "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.bam",
        bai = "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.bam.bai"
    output:
        bam = "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.mkdup.bam",
        metrics = "stats/{female_sample}.GRCh38.p12.genome.XXonly.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule index_mkdup_bam_females:
    input:
        "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.mkdup.bam"
    output:
        "processed_bams/{female_sample}.GRCh38.p12.genome.XXonly.sorted.mkdup.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"
