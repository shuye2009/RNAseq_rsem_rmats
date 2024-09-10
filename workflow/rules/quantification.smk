resultdir = config['path']['resultdir']

rule build_transcriptome:
    input:
        genome = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        'resources/transcriptome.fa'
    log:
        resultdir+"/logs/gffread.log"
    message:
        "Build a reference transcriptome using gffread."
    conda:
        "confflinks-2.2.1"
    resources:
        mem_mb = 4000,
        runtime = 60
    shell:
        "gffread {input.gtf} -g {input.genome} -w {output} &> {log}"


rule kallisto_index:
    input:
        rules.build_transcriptome.output
    output:
        resultdir+"/kallisto/kallisto.idx"
    params:
        31
    log:
        resultdir+"/logs/kallisto/index.log"
    message:
        "Build a Kallisto index from the transcriptome FASTA file."
    conda:
        "kallisto-0.48.0"
    resources:
        mem_mb = 4000,
        runtime = 120
    shell:
        "kallisto index "
        "--index={output} "
        "--kmer-size={params} "
        "{input} "
        "&> {log}"


rule kallisto_quant:
    input:
        idx = rules.kallisto_index.output,
        fq1 = rules.trimmomatic.output.r1,
        fq2 = rules.trimmomatic.output.r2
    output:
        tsv = resultdir+"/kallisto/{sample}/abundance.tsv",
        h5 = resultdir+"/kallisto/{sample}/abundance.h5"
    params:
        bootstrap = "50",
        outdir = resultdir+"/kallisto/{sample}"
    threads:
        10
    log:
        resultdir+"/logs/kallisto/{sample}.log"
    message:
        "Perform pseudoalignment and quantify transcript abundance for {wildcards.sample}."
    conda:
        "kallisto-0.48.0"
    resources:
        mem_mb = 16000,
        runtime = 120
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"


rule tx2gene:
    input:
        gtf = config["path"]["genome_gtf"]
    output:
        tsv = "resources/tx2gene.tsv",
        gene_tsv = "resources/geneid2name.tsv"
    message:
        "Convert transcript IDs to gene IDs."
    conda:
        "polars-0.20.7"
    resources:
        mem_mb = 3000,
        runtime = 60
    script:
        "../scripts/tx2gene.py"
       

## Isoform level merged table
rule merge_kallisto_quant:
    input:
        quant = expand(rules.kallisto_quant.output.tsv, sample=SAMPLES),
        tx2gene = rules.tx2gene.output.tsv,
        gtf = config["path"]["genome_gtf"]
    output:
        tpm = resultdir+"/kallisto/tpm_kallisto.tsv"
    message: 
        "Merge kallisto quantification results into one dataframe for further analysis."
    conda:
        "polars-0.20.7"
    resources:
        mem_mb = 8000,
        runtime = 120
    script:
        "../scripts/merge_kallisto_quant.py"