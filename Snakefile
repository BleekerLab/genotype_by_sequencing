#################################################
# Snakemake pipeline for genotyping by sequencing 
#################################################


###########
# Libraries
###########
import pandas as pd
import subprocess

###############
# Configuration
###############
configfile: "config/config.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]


########################
# Samples and conditions
########################

# create lists containing the sample names and conditions
samples = pd.read_csv(config["samples"], dtype=str,index_col=0,sep="\t")
SAMPLES = samples.index.get_level_values('sample').unique().tolist()


###########################
# Input functions for rules
###########################

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the samples.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """This function checks if the sample has paired end or single end reads and returns 1 or 2 names of the fastq files"""
    if sample_is_single_end(wildcards.sample):
        return samples.loc[(wildcards.sample), ["fq1"]].dropna()
    else:
        return samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_trim_names(wildcards):
    """
    This function:
      1. Checks if the sample is paired end or single end
      2. Returns the correct input and output trimmed file names. 
    """
    if sample_is_single_end(wildcards.sample):
        inFile = samples.loc[(wildcards.sample), ["fq1"]].dropna()
        return "--in1 " + inFile[0] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq" 
    else:
        inFile = samples.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
        return "--in1 " + inFile[0] + " --in2 " + inFile[1] + " --out1 " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq --out2 "  + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq"

def get_star_names(wildcards):
    """
    This function:
      1. Checks if the sample is paired end or single end.
      2. Returns the correct input file names for STAR mapping step.
    """
    if sample_is_single_end(wildcards.sample):
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz"     
    else:
        return WORKING_DIR + "trimmed/" + wildcards.sample + "_R1_trimmed.fq.gz " + WORKING_DIR + "trimmed/" + wildcards.sample + "_R2_trimmed.fq.gz"

#################
# Desired outputs
#################
MULTIQC = RESULT_DIR + "multiqc_report.html"
MAPPING_REPORT = RESULT_DIR + "mapping_summary.csv"
SNP_COUNTS = expand(RESULT_DIR + "{sample}.counts.tsv", sample = SAMPLES)

if config["keep_working_dir"] == True:
    rule all:
        input:
            MULTIQC,
            MAPPING_REPORT,
            SNP_COUNTS   
        message:
            "Genotyping by sequencing pipeline run complete!"
        shell:
            "cp config/config.yaml {RESULT_DIR};"
            "cp config/samples.tsv {RESULT_DIR};"
else:
    rule all:
        input:
            MULTIQC,
            MAPPING_REPORT,
            SNP_COUNTS   
        message:
            "Genotyping by sequencing pipeline run complete!"
        shell:
            "cp config/config.yaml {RESULT_DIR};"
            "cp config/samples.tsv {RESULT_DIR};"
            "rm -r {WORKING_DIR}"

#######
# Rules
#######


###########################
# Genome reference indexing
###########################

if config["datatype"] == "RNA":
    rule star_index:
        input:
            fasta = config["refs"]["genome"],
            gtf =   config["refs"]["gtf"]
        output:
            genome_index = [WORKING_DIR + "genome/" + f for f in ["chrLength.txt","chrNameLength.txt","chrName.txt","chrStart.txt","Genome","genomeParameters.txt","SA","SAindex"]]
        message:
            "generating STAR genome index"
        params:
            genome_dir = WORKING_DIR + "genome/",
            sjdb_overhang = config["star_index"]["sjdbOverhang"],
            limit_genome_generate_ram = config["star_index"]["limitGenomeGenerateRAM"],
            genome_sa = config["star_index"]["genomeSAindexNbases"],
            genome_chr_bin_n_bits = config["star_index"]["genomeChrBinNbits"]
        threads: 20
        resources: mem_mb=100000
        shell:
            "mkdir -p {params.genome_dir}; " # if directory not created STAR will ask for it
            "STAR --runThreadN {threads} "
            "--runMode genomeGenerate "
            "--genomeDir {params.genome_dir} "
            "--genomeFastaFiles {input.fasta} "
            "--sjdbGTFfile {input.gtf} "
            "--sjdbOverhang {params.sjdb_overhang} "
            "--limitGenomeGenerateRAM {params.limit_genome_generate_ram} "
            "--genomeSAindexNbases {params.genome_sa} "
            "--genomeChrBinNbits {params.genome_chr_bin_n_bits}"
elif config["datatype"] == "DNA":
    rule bwa_index:
        input:
            fasta = config["refs"]["genome"]
        output: 
            genome_index = [WORKING_DIR + "genome/genome" + ext for ext in [".sa", ".pac", ".bwt", ".ann", ".amb"] ]
        message:
            "Generating BWA genome index"
        params:
            genome_dir = WORKING_DIR + "genome/"
        shell:
            "bwa index -p genome {input.fasta};"
            "mv genome.amb genome.ann genome.bwt genome.pac genome.sa {params.genome_dir}"
else:
    raise ValueError('Please specify either "DNA" or "RNA" as "datatype" in the config.yaml file.')


#######################
# RNA-seq read trimming
#######################

rule fastp:
    input:
        get_fastq
    output:
        fq1  = temp(WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq"),
        fq2  = temp(WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq"),
        html = WORKING_DIR + "fastp/{sample}_fastp.html",
        json = WORKING_DIR + "fastp/{sample}_fastp.json"
    message:"trimming {wildcards.sample} reads"
    threads: 20
    log:
        RESULT_DIR + "fastp/{sample}.log.txt"
    params:
        sampleName = "{sample}",
        in_and_out_files =  get_trim_names,
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    resources: cpus=10
    shell:
        "touch {output.fq2};\
        fastp --thread {threads}  --html {output.html} --json {output.json} \
        --qualified_quality_phred {params.qualified_quality_phred} \
        {params.in_and_out_files} \
        2>{log}"

rule multiqc:
    input:
        expand(WORKING_DIR + "fastp/{sample}_fastp.json", sample = SAMPLES)
    output:
        RESULT_DIR + "multiqc_report.html"
    params:
        fastp_directory = WORKING_DIR + "fastp/",
        outdir = RESULT_DIR
    message: "Summarising fastp reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "

####################################
# DNA-seq or RNA-Seq read alignement
####################################

if config["datatype"] == "RNA":
    rule map_to_genome_using_STAR:
        input:
            genome_index = rules.star_index.output,
            forward_read = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq",
            reverse_read = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq"
        output:
            WORKING_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam",
            WORKING_DIR + "star/{sample}_Log.final.out"
        message:
            "mapping {wildcards.sample} RNA-seq reads to genome"
        params:
            sample_name           =  "{sample}",
            star_input_file_names =  get_star_names,
            prefix                =  WORKING_DIR + "star/{sample}_",
            maxmismatches         =  config["star"]["mismatches"],
            unmapped              =  config["star"]["unmapped"]   ,
            multimappers          =  config["star"]["multimappers"],
            matchNminoverLread    =  config["star"]["matchminoverlengthread"],
            outSamType            =  config["star"]["samtype"],
            outSAMattributes      =  config["star"]["samattributes"],
            intronmax             =  config["star"]["intronmax"],
            matesgap              =  config["star"]["matesgap"],
            genome_index          =  WORKING_DIR + "genome/"
        threads: 20
        shell:
            """
            STAR --genomeDir {params.genome_index} --readFilesIn {params.star_input_file_names} --readFilesCommand zcat --outFilterMultimapNmax {params.multimappers} \
            --outFilterMismatchNmax {params.maxmismatches} --alignMatesGapMax {params.matesgap} --alignIntronMax {params.intronmax}  \
            --outFilterMatchNminOverLread  {params.matchNminoverLread} --alignEndsType EndToEnd --runThreadN {threads}  --outReadsUnmapped {params.unmapped} \
            --outFileNamePrefix {params.prefix} --outSAMtype {params.outSamType}  --outSAMattributes {params.outSAMattributes}
            """
elif config["datatype"] == "DNA":
    rule map_to_genome_using_bwa:
        input:
            genome_index = rules.bwa_index.output,
            forward_read = WORKING_DIR + "trimmed/" + "{sample}_R1_trimmed.fq",
            reverse_read = WORKING_DIR + "trimmed/" + "{sample}_R2_trimmed.fq"
        output:
            WORKING_DIR + "bwa/{sample}_aligned.sorted.bam"
        message: 
            "Mapping {wildcards.sample} DNA-seq reads to genome"
        threads: 20
        params:
           genome_index = WORKING_DIR + "genome/genome"
        shell:
            "bwa mem -v 0 -t {threads} {params.genome_index} {input.forward_read} {input.reverse_read} | samtools sort -@ {threads} -o {output} - "
          
if config["datatype"] == "RNA":
    rule generate_mapping_summary:
        input:
            expand(WORKING_DIR + "star/{sample}_Log.final.out", sample = SAMPLES)
        output:
            RESULT_DIR + "mapping_summary.csv"
        message:
            "Concatenating STAR mapping report from RNA-seq data and generating .csv mapping summary."
        params:
            directory_with_mapping_reports = WORKING_DIR + "star/",
            config_file_path = "config/config.yaml"
        shell:
            "python scripts/generate_mapping_summary.py {params.directory_with_mapping_reports} {params.config_file_path} {output}"
elif config["datatype"] == "DNA":
    rule generate_mapping_summary:
        input:
            expand(WORKING_DIR + "bwa/{sample}_aligned.sorted.bam", sample = SAMPLES)
        output: 
            RESULT_DIR + "mapping_summary.csv"
        message:
            "Creating BWA mapping report from DNA-seq data and generate .csv mapping summary"
        params: 
            directory_with_bam_files = WORKING_DIR + "bwa/"
        shell:
            "touch {output}"


########################
# SNP calling per sample
########################

if config["datatype"] == "RNA":
    rule call_snps:
        input:
            bam = WORKING_DIR + "star/{sample}_Aligned.sortedByCoord.out.bam",
            fasta = config["refs"]["genome"],
        output:
            WORKING_DIR + "vcf/{sample}.vcf"
        message:
            "Calling SNPs using bcftools from {wildcards.sample} mapping"
        params: 
            max_depth = config["bcftools"]["max_depth"],
            min_base_quality = config["bcftools"]["min_base_quality"]
        threads: 20
        shell:
            "bcftools mpileup --threads {threads} "
            "--max-depth {params.max_depth} "
            "--no-BAQ " # Disable probabilistic realignment for the computation of base alignment quality (BAQ). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments.
            "--min-BQ {params.min_base_quality} "
            "--per-sample-mF "
            "-f {input.fasta} "
            "--skip-indels "
            "--output-type v {input.bam} | "
            "bcftools call --multiallelic-caller --variants-only -Ov - > {output}" 

if config["datatype"] == "DNA":
    rule call_snps:
        input:
            bam = WORKING_DIR + "bwa/{sample}_aligned.sorted.bam",
            fasta = config["refs"]["genome"],
        output:
            WORKING_DIR + "vcf/{sample}.vcf"
        message:
            "Calling SNPs using bcftools from {wildcards.sample} mapping"
        params: 
            max_depth = config["bcftools"]["max_depth"],
            min_base_quality = config["bcftools"]["min_base_quality"]
        threads: 20
        shell:
            "bcftools mpileup --threads {threads} "
            "--max-depth {params.max_depth} "
            "--no-BAQ " # Disable probabilistic realignment for the computation of base alignment quality (BAQ). BAQ is the Phred-scaled probability of a read base being misaligned. Applying this option greatly helps to reduce false SNPs caused by misalignments.
            "--min-BQ {params.min_base_quality} "
            "--per-sample-mF "
            "-f {input.fasta} "
            "--skip-indels "
            "--output-type v {input.bam} | "
            "bcftools call --multiallelic-caller --variants-only -Ov - > {output}" 

rule filter_snps_based_on_quality:
    input: 
        WORKING_DIR + "vcf/{sample}.vcf"
    output:
        WORKING_DIR + "vcf/{sample}.qual.vcf"
    message:
        "Filtering {wildcards.sample} SNPs based on quality"
    params:
        snp_quality = config["bcftools"]["snp_quality"],
        read_depth = config["bcftools"]["read_depth"]
    threads: 20
    shell:
        "bcftools view "
        "--output-type v "
        "--output-file {output} "
        "--threads {threads} "
        "--include 'MIN(DP)>{params.read_depth} && MIN(QUAL)>{params.snp_quality}' "
        "{input}"


rule keep_only_homozygous_alt_genotypes:
    input:
        WORKING_DIR + "vcf/{sample}.qual.vcf"
    output:
        WORKING_DIR + "vcf/{sample}.qual.alt.vcf"
    message:
        "Keeping homozygous ALT-ALT homozygous genotypes from {wildcards.sample} VCF file"
    threads: 20
    params:
        bcftools_filter_expr = "GT='AA'"
    shell:
        "bcftools view --output-type v "
        "-o {output} "
        "--include {params.bcftools_filter_expr:q} " # robust quoting https://carpentries-incubator.github.io/snakemake-novice-bioinformatics/13-quoting/index.html
        "{input}"

rule convert_vcf_to_bed:
    input:
        WORKING_DIR + "vcf/{sample}.qual.alt.vcf"
    output:
        WORKING_DIR + "bed/{sample}.bed"
    message:
        "Convert {wildcards.sample} VCF to BED format"
    shell:
        "vcf2bed < {input} > {output}"


####################
# Create genome bins
####################

rule compute_chromosome_sizes:
    input:
        fasta = config["refs"]["genome"]
    output:
        WORKING_DIR + "genome/chromsizes.txt"
    message:
        "Compute genome chromosome sizes"
    params: 
        samtools_index_file_name = config["refs"]["genome"] + ".fai"
    shell:
        "samtools faidx {input.fasta} ; cut -f1,2 {params.samtools_index_file_name} > {output}" 

rule parse_chromosome_sizes:
    input: 
        WORKING_DIR + "genome/chromsizes.txt"
    output:
        WORKING_DIR + "genome/chromsizes.parsed.txt"
    message:
        "Parse chromosome sizes file to keep the two first columns"
    shell:
        "cut -f 1,2 {input} > {output}"

rule create_genome_bins:
    input:
        chromsizes = WORKING_DIR + "genome/chromsizes.parsed.txt"
    output:
        bed = WORKING_DIR + "genome/genome.bed"
    message:
        "Create a BED file from genome chromosome sizes of size {params.window_size}"
    params:
        window_size = config["bedtools"]["window_size"]
    shell:
        "bedtools makewindows -g {input} -w {params.window_size} > {output}"

##################################
# Create SNP counts per genome bin
##################################

rule count_nb_snp_per_genome_bin:
    input:
        sample_bed = WORKING_DIR + "bed/{sample}.bed",
        genome_bed = WORKING_DIR + "genome/genome.bed"
    output:
        RESULT_DIR + "{sample}.counts.tsv"
    message:
        "Count number of SNPs in {wildcards.sample} per {params.window_size}"
    params:
        window_size = config["bedtools"]["window_size"] 
    shell:
        "bedmap --echo --count --delim '\t' {input.genome_bed} {input.sample_bed} > {output}" 

