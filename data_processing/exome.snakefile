import os

configfile: "exome_config.json"

# Tool paths:
fastqc_path = "fastqc"
multiqc_path = "multiqc"
bwa_path = "bwa"
samtools_path = "samtools"
bbduksh_path = "bbduk.sh"

#
REF_TYPE = ["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]

# Directory
fastq_directory = "fastq_files/"

rule all:
    input:
    	"multiqc_results/multiqc_report.html"
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

rule fastqc_analysis:
	input:
		"fastq_files/{sample_name}_{read}.fastq.gz"
	output:
		"fastqc_results/{sample_name}.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"PERL5LIB=/home/tphung3/softwares/miniconda3/envs/epitopepipeline/lib/site_perl/5.26.2/ {params.fastqc} -o fastqc_results {input}"

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
