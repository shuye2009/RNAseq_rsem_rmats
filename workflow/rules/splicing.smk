resultdir = config['path']['resultdir']


rule make_files:
    input:
        files = expand(rules.primary_alignments.output, sample=SAMPLES),
        design = "config/design.tsv",
        comps = "config/comparisons.tsv"
    output:
        b1 = expand("config/{comp}/rmats_b1.txt", comp=comparisons),
        b2 = expand("config/{comp}/rmats_b2.txt", comp=comparisons)
    message:
        "Create b1 b2 files for rmats."
    script:
        "../scripts/group_files.py"

rule rmats:
    input:
        bam = expand(rules.primary_alignments.output, sample=SAMPLES),
        group1 = "config/{comp}/rmats_b1.txt",
        group2 = "config/{comp}/rmats_b2.txt",
        gtf = config["path"]["genome_gtf"]
    output:
        outdir = directory(resultdir+"/rmats/{comp}/raw"),
        tmpdir = directory(resultdir+"/rmats/{comp}/tmp"),
        summary = resultdir+"/rmats/{comp}/raw/summary.txt"
    params:
        readlength = config["params"]["readlength"]
    threads:
        4
    message:
        "Run rMATS for {wildcards.comp}."
    conda:
        "rmats-4.3.0"
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        "rmats.py --b1 {input.group1} --b2 {input.group2} "
        "--gtf {input.gtf} -t paired --readLength {params.readlength} --variable-read-length "
        "--nthread {threads} --od {output.outdir} --tmp {output.tmpdir} "
        "--anchorLength 1 --task both --novelSS"


rule filter_rmats:
    input:
        summary = rules.rmats.output.summary,
        tpm = rules.merge_rsem.output.tpm
    output:
        result = resultdir+"/rmats/{comp}/filtered/SE.tsv"
    params:
        dir = directory(resultdir+"/rmats/{comp}"),
        gtf = config["path"]["genome_gtf"],
        fdr = 0.05,
        deltapsi = 0.10
    message:
        "Filter raw rMATS output for {wildcards.comp}."
    conda:
        "polars-0.20.7"
    resources:
        mem_mb = 4000,
        runtime = 60
    script:
        "../scripts/filter_rmats.py"
