sampledir = config['path']['datadir']
resultdir = config['path']['resultdir']

rule fastqc_pretrim:
    input:
        sampledir+"/{sample}_{read}.fastq.gz"
    output:
        resultdir+"/fastqc/pretrim/{sample}_{read}_fastqc.html"
    params:
        resultdir+"/fastqc/pretrim"
    log:
        resultdir+"/logs/fastqc/pretrim/{sample}_{read}.log"
    threads:
        32
    message:
        "Quality control check on raw sequence data of {wildcards.sample}_{wildcards.read}."
    conda:
        "fastqc-0.12.1"
    resources:
        mem_mb = 8000,
        runtime = 60,
        partition = "himem"
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "-t {threads} "
        "{input} "
        "&> {log}"


rule trimmomatic:
    input:
        r1 = sampledir+"/{sample}_R1.fastq.gz",
        r2 = sampledir+"/{sample}_R2.fastq.gz",
    output:
        r1 = resultdir+"/trimmomatic/trimmed_{sample}_R1.fastq.gz",
        r2 = resultdir+"/trimmomatic/trimmed_{sample}_R2.fastq.gz",
        unpaired_r1 = resultdir+"/trimmomatic/{sample}_R1.unpaired.fastq.gz",
        unpaired_r2 = resultdir+"/trimmomatic/{sample}_R2.unpaired.fastq.gz"
    log:
        resultdir+"/logs/trimmomatic/{sample}.log"
    params:
        trimmer = "ILLUMINACLIP:resources/Adapters-PE_NextSeq.fa:2:12:10:8:true TRAILING:30 LEADING:30 MINLEN:20",
    threads:
        8
    message:
        "Trim poor quality reads in {wildcards.sample} using Trimmomatic."
    conda:
        "trimmomatic-0.39"
    resources:
        mem_mb = 16000,
        runtime = 60,
        partition = "himem"
    shell:
        "trimmomatic PE "
        "-threads {threads} "
        "-phred33 "
        "{input.r1} {input.r2} "
        "{output.r1} {output.unpaired_r1} "
        "{output.r2} {output.unpaired_r2} "
        "{params.trimmer} "
        "2> {log}"


rule fastqc_posttrim:
    input:
        r = resultdir+"/trimmomatic/trimmed_{sample}_{read}.fastq.gz"
    output:
        html = resultdir+"/fastqc/posttrim/trimmed_{sample}_{read}_fastqc.html"
    params:
        resultdir+"/fastqc/posttrim"
    threads:
        32
    log:
        resultdir+"/logs/fastqc/posttrim/trimmed_{sample}_{read}.log"
    message:
        "Quality control check on trimmed sequence data of {wildcards.sample}_{wildcards.read}."
    conda:
        "fastqc-0.12.1"
    resources:
        mem_mb = 8000,
        runtime = 60,
        partition = "himem"
    shell:
        "fastqc "
        "--outdir {params} "
        "--format fastq "
        "-t {threads} "
        "{input.r} "
        "&> {log}"
    