sampledir = config['path']['datadir']
resultdir = config['path']['resultdir']

rule star_index:
    input:
        fasta = config["path"]["genome_fasta"],
        gtf = config["path"]["genome_gtf"]
    output:
        config["path"]["indexdir"]+"/chrNameLength.txt"
    params:
        idx = config["path"]["indexdir"]
    threads: 
        32
    log:
        resultdir+"/logs/STAR/index.log"
    message:
        "Generate genome indexes files using STAR."
    conda:
        "star-2.7.11b"
    resources:
        partition = "himem",
        mem_mb = 32000,
        runtime = 600

    shell:
        "STAR --runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir {params.idx} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang 99 "
        "&> {log}"



rule star_align:
    input:
        fq1 = rules.trimmomatic.output.r1,
        fq2 = rules.trimmomatic.output.r2,
        chrNameLength = rules.star_index.output
    output:
        bam = resultdir+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    log:
        resultdir+"/logs/STAR/{sample}_align.log"
    params:
        out_prefix = resultdir+"/STAR/{sample}/{sample}_",
        idx = rules.star_index.params.idx
    threads:
        32
    message:
        "Align {wildcards.sample} reads to the reference genome using STAR."
    conda:
        "star-2.7.11b"
    resources:
        partition = "himem",
        mem_mb = 32000,
        runtime = 600
    shell:
        "STAR --runMode alignReads "
        "--genomeDir {params.idx} "
        "--readFilesIn {input.fq1} {input.fq2} "
        "--runThreadN {threads} "
        "--readFilesCommand zcat --outReadsUnmapped Fastx "
        "--outFilterType BySJout --outStd Log --outSAMunmapped Within "
        "--outSAMtype BAM SortedByCoordinate "
        "--outFileNamePrefix {params.out_prefix} "
        "--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 "
        "--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 "
        "--outSAMstrandField intronMotif "
        "--outSAMattributes NH HI AS NM MD XS --twopassMode Basic "
        "--alignEndsProtrude 5 ConcordantPair "
        "--alignEndsType EndToEnd " 
        "--alignSJDBoverhangMin 1 --alignSJoverhangMin 8 "
        "&> {log}"

rule primary_alignments:
    input:
        rules.star_align.output.bam
    output:
        resultdir+"/STAR/{sample}/{sample}_Aligned.sortedByCoord.out.primary.bam"
    log:
        resultdir+"/logs/STAR/{sample}_primary.log"
    message:
        "Keep primary alignments only for {wildcards.sample}."
    conda:
        "star-2.7.11b"
    resources:
        mem_mb = 3000,
        runtime = 60
    shell:
        "samtools view -b -F 256 -o {output} {input} "
        "&> {log}"

rule picard:
    input:
        rules.primary_alignments.output
    output:
        txt = resultdir+"/picard/{sample}.isize.txt",
        pdf = resultdir+"/picard/{sample}.isize.pdf"
    log:
        resultdir+"/logs/picard/{sample}.log"
    message:
        "Collect insert size distribution metrics for validating library construction for {wildcards.sample}."
    conda:
        "picard-3.1.1"
    resources:
        mem_mb = 8000,
        runtime = 60
    shell:
        "picard CollectInsertSizeMetrics "
        "--INPUT {input} "
        "--OUTPUT {output.txt} "
        "--Histogram_FILE {output.pdf} "
        "--MINIMUM_PCT 0.5 "
        "&> {log}"


rule genomecov:
    input:
        rules.primary_alignments.output
    output:
        resultdir+"/genomecov/{sample}.bedgraph"
    message:
        "Report {wildcards.sample} genome coverage in BEDGRAPH format."
    conda:
        "bedtools-2.31.1"
    resources:
        mem_mb = 8000,
        runtime = 60
    shell:
        "bedtools genomecov -bg -split -ibam {input} | sort -k1,1 -k2,2n > {output}"