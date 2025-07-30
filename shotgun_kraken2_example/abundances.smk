# Snakemake workflow to use kraken2 + bracken to calculate species-level
# abundances in the provided samples of paired-end reads.

# Invoke snakemake with --config datadir=<path_to_data>
# If not provided, input data files are assumed to exist in the current working # directory.
try:
    exp_datadir = config["datadir"]
    datadir = exp_datadir
except(KeyError):
    exp_datadir = "."
    datadir = ""

SAMPLES = [d for d in os.listdir(exp_datadir) if d.startswith('ERR')]

rule all:
    input:
        expand("abund/bracken/{sample}.bracken", sample=SAMPLES)

rule kraken:
    input: 
        in1 = os.path.join(datadir, "{sample}/{sample}.sra_1.fastq"),
        in2 = os.path.join(datadir, "{sample}/{sample}.sra_2.fastq"),
    output: 
        out = "abund/{sample}.kraken2",
        report = "abund/{sample}.kreport",
    resources:
        time = "0:10:00",
        mem_mb = "100000",
        nodes = 1,
    shell:
        "kraken2 --db /fs/ess/PAS2276/db/uhgg_kraken2-db "
        "--threads 40 "
        "--output {output.out} "
        "--report {output.report} "
        "--paired "
        "{input.in1} {input.in2} "


rule bracken:
    input: 
        "abund/{sample}.kreport",
    output: 
        "abund/{sample}.bracken",
    resources:
        time = "0:10:00",
        nodes = 1,
    shell:
        "bracken -d /fs/ess/PAS2276/db/uhgg_kraken2-db "
        "-i {input} "
        "-o {output} "
        "-r 100 "
        "-l R7"
