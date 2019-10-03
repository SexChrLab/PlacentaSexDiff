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

# ----------------------------------
# Exome processing for males samples
# ----------------------------------
rule all:
    input:
        expand("processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.mkdup.bam.bai", male_sample = config["males"])

rule prep_refs_mk_sy_ln:
    input:
        fa = config["XY_ref_fa"],
        fai = config["XY_ref_fai"]
    output:
        fa = "refs/GRCh38.p12.genome.XY.fa",
        fai = "refs/GRCh38.p12.genome.XY.fa.fai"
    shell:
        """
        ln -s {input.fa} {output.fa};
        ln -s {input.fai} {output.fai}
        """

rule prep_refs_males:
    input:
        ref = "refs/GRCh38.p12.genome.XY.fa"
    output:
        amb = "refs/GRCh38.p12.genome.XY.fa.amb",
        dict = "refs/GRCh38.p12.genome.XY.dict"
    params:
        samtools = samtools_path,
        bwa = bwa_path
    run:
        # .dict
        shell("{params.samtools} dict -o {output.dict} {input.ref}")

        # bwa
        shell("{params.bwa} index {input.ref}")

rule map_to_sex_specific_refs_males:
    input:
        fq1 = "/data/storage/SAYRES/placenta_YPOPS/08_exome/trimmed_fastqs/{male_sample}_trimmed_R1.fastq.gz",
        fq2 = "/data/storage/SAYRES/placenta_YPOPS/08_exome/trimmed_fastqs/{male_sample}_trimmed_R2.fastq.gz",
        ref = "refs/GRCh38.p12.genome.XY.fa",
        fai = "refs/GRCh38.p12.genome.XY.fa.fai"
    output:
        "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.bam"
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

rule index_bam_males:
    input:
        "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.bam"
    output:
        "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule picard_mkdups_males:
    input:
        bam = "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.bam",
        bai = "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.bam.bai"
    output:
        bam = "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.mkdup.bam",
        metrics = "stats/{male_sample}.GRCh38.p12.genome.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
        "M={output.metrics}"

rule index_mkdup_bam_males:
    input:
        "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.mkdup.bam"
    output:
        "processed_bams/{male_sample}.GRCh38.p12.genome.XY.sorted.mkdup.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"
