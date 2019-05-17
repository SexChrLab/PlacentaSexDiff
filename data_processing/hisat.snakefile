import os

configfile: "RNA_config.json"

# Tool paths:
fastqc_path = "fastqc"
multiqc_path = "multiqc"
bwa_path = "bwa"
samtools_path = "samtools"
bbduksh_path = "bbduk.sh"
hisat_path = "hisat2"
bamtools_path = "bamtools"
picard_path = "picard"
featureCounts_path = "featureCounts"


#
REF_TYPE = ["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]
REF_TYPE_HISAT = ["Ref_GRCh38_Y_HardMasked_HISAT_index","Ref_GRCh38_Y_PARsMasked_HISAT_index"]

# Directory
fastq_directory = "fastq_files/"

rule all:
    input:
        expand(config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_XX.sam", female_sample = config["females"]),
        expand(config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_XY.sam", male_sample = config["males"]),
        
        expand(config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_XX.bam", female_sample = config["females"]),
        expand(config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_XY.bam", male_sample = config["males"]),
        
        expand(config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam", male_sample = config["males"]),
        expand(config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam", female_sample = config["females"]),
        
        expand(config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam", male_sample = config["males"]),
        expand(config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam", female_sample = config["females"]),
        
        expand(config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam", male_sample = config["males"]),
        expand(config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam", female_sample = config["females"]),

        expand(config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam.bai", male_sample = config["males"]),
        expand(config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam.bai", female_sample = config["females"]),
        
        #expand(config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XY.txt", male_sample = config["males"]),
        #expand(config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XX.txt", female_sample = config["females"]),

        expand("FeatureCounts/{male_sample}_HISAT_geneCounts_XY.txt", male_sample = config["males"]),
        expand("FeatureCounts/{female_sample}_HISAT_geneCounts_XX.txt", female_sample = config["females"]),

        expand("FeatureCounts/{male_sample}_HISAT_transcriptCounts_XY.txt", male_sample = config["males"]),
        expand("FeatureCounts/{female_sample}_HISAT_transcriptCounts_XX.txt", female_sample = config["females"])

    input:
        expand("fastqc_results/{sample_name}_trimmed_R1.fastqc.html", sample_name =config["sample_names"]),
        expand("fastqc_results/{sample_name}_trimmed_R2.fastqc.html", sample_name =config["sample_names"])
    input:
        expand("trimmed_fastqs/{sample_name}_trimmed_R1.fastq", sample_name = config["sample_names"]),
        expand("trimmed_fastqs/{sample_name}_trimmed_R2.fastq", sample_name = config["sample_names"])
    input:
        expand("fastq_files/{sample_name}_R1.fastq", sample_name = config["sample_names"])
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
        R1_out = "fastq_files/{sample_name}_R1.fastq",
        R2_out = "fastq_files/{sample_name}_R2.fastq"
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
        out_fq1 = "trimmed_fastqs/{sample_name}_trimmed_R1.fastq",
        out_fq2 = "trimmed_fastqs/{sample_name}_trimmed_R2.fastq"
    params:
        bbduksh = bbduksh_path
    shell:
        "{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
        "out1={output.out_fq1} out2={output.out_fq2} "
        "ref=/mnt/storage/SAYRES/REFERENCE_GENOMES/adapters/adapter_sequence.fa "
        "qtrim=rl trimq=30 minlen=75 maq=20"

rule HISAT_paired_males:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{male_sample}_trimmed_R1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{male_sample}_trimmed_R2.fastq"
    output:
        out_1 = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_XY.sam")
    params:
        HISAT_Index_male = config["HG38_Transcriptome_Index_HISAT_Path_male"],
    shell:
        "hisat2 --dta -q --phred33 -p 8 -x {params.HISAT_Index_male} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"

rule HISAT_paired_females:
    input:
        Trimmed_FASTQ1 = "trimmed_fastqs/{female_sample}_trimmed_R1.fastq",
        Trimmed_FASTQ2 = "trimmed_fastqs/{female_sample}_trimmed_R2.fastq"
    output:
        out_1 = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_XX.sam")
    params:
        HISAT_Index_female = config["HG38_Transcriptome_Index_HISAT_Path_female"],
    shell:
        "hisat2 --dta -q --phred33 -p 8 -x {params.HISAT_Index_female} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"


# males
rule samtools_view_males:
    input:
        SAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_XY.sam")
    output:
        BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_XY.bam")
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule bam_sort_males:    
    input:
        IN_BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_XY.bam")
    output:
        sort_BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_XY.bam")
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups_males:
    input:
        sort_BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_XY.bam")
    output:
        BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam"),
        metrics = "stats/{male_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps_males:
    input:
        Read_BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam")
    output:
        BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"),
    params:
        id = lambda wildcards: config[wildcards.male_sample]["ID"],
        sm = lambda wildcards: config[wildcards.male_sample]["SM"],
        lb = lambda wildcards: config[wildcards.male_sample]["LB"],
        pu = lambda wildcards: config[wildcards.male_sample]["PU"],
        pl = lambda wildcards: config[wildcards.male_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam_males:
    input:
        BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam")
    output:
        BAI = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam.bai")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule stats_bam_males:
    input:
        BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam")
    output:
        stats = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XY.txt")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule feautreCounts_gene_males:
    input:
        BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam")
    output:
        counts = "FeatureCounts/{male_sample}_HISAT_geneCounts_XY.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 0 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

rule feautreCounts_transcripts_males:
    input:
        BAM = (config["Output_File_Directory"]+"{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam")
    output:
        counts = "FeatureCounts/{male_sample}_HISAT_transcriptCounts_XY.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 0 -t exon -g transcript_id -O -a {params.GTF} -o {output.counts} {input.BAM}"
        
# females
rule samtools_view_females:
    input:
        SAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_XX.sam")
    output:
        BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_XX.bam")
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule bam_sort_females:  
    input:
        IN_BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_XX.bam")
    output:
        sort_BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_XX.bam")
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups_females:
    input:
        sort_BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_XX.bam")
    output:
        BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam"),
        metrics = "stats/{female_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps_females:
    input:
        Read_BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam")
    output:
        BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"),
    params:
        id = lambda wildcards: config[wildcards.female_sample]["ID"],
        sm = lambda wildcards: config[wildcards.female_sample]["SM"],
        lb = lambda wildcards: config[wildcards.female_sample]["LB"],
        pu = lambda wildcards: config[wildcards.female_sample]["PU"],
        pl = lambda wildcards: config[wildcards.female_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam_females:
    input:
        BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam")
    output:
        BAI = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam.bai")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule stats_bam_females:
    input:
        BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam")
    output:
        stats = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XX.txt")
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"
 
rule feautreCounts_gene_females:
    input:
        BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam")
    output:
        counts = "FeatureCounts/{female_sample}_HISAT_geneCounts_XX.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 0 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"
        
rule feautreCounts_transcripts_females:
    input:
        BAM = (config["Output_File_Directory"]+"{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam")
    output:
        counts = "FeatureCounts/{female_sample}_HISAT_transcriptCounts_XX.txt"
    params:
        featureCounts = featureCounts_path,
        GTF = config["gencode.v29.annotation.gtf_path"],
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 0 -t exon -g transcript_id -O -a {params.GTF} -o {output.counts} {input.BAM}"
