### Adapted from Danny Bergeron's code
resultdir = config['path']['resultdir']

rule gene_abundance:
    input:
        quant = expand(rules.kallisto_quant.output, sample=SAMPLES),
        samples = "config/design.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        out_file = resultdir+"/deseq2/gene_level_abundance.tsv"
    params:
        kallisto_dir = resultdir+"/kallisto",
        out_dir = resultdir+"/deseq2"
    message:
        "Perform gene level abundance calculation."
    conda:
        "deseq2-1.42.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/gene_level_abundance.R"

rule deseq2:
    input:
        quant = expand(rules.kallisto_quant.output, sample=SAMPLES),
        samples = "config/design.tsv",
        comparisons = "config/comparisons.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        out_files = resultdir+"/deseq2/{comp}.csv"
    params:
        kallisto_dir = resultdir+"/kallisto",
        out_dir = resultdir+"/deseq2"
    message:
        "Perform differential expression analysis for {wildcards.comp}."
    conda:
        "deseq2-1.42.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/DESeq2_kallisto_tximport.R"


rule volcano_plot:
    input:
        DE_output = rules.deseq2.output.out_files,
        filtered_genes = rules.merge_kallisto_quant.output.tpm,
    output:
        volcano = resultdir+"/deseq2/{comp}.svg",
        up_genes = resultdir+"/deseq2/{comp}_sig_DE_up.tsv",
        down_genes = resultdir+"/deseq2/{comp}_sig_DE_down.tsv"
    params:
        pval_threshold = 0.05,
        gtf = config["path"]["genome_gtf"]
    message:
        "Create a volcano plot using deseq2 output for {wildcards.comp}."
    conda:
        "gprofiler-1.0.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/volcano_plot.py"


rule GO_upregulated_genes:
    input:
        genes = rules.volcano_plot.output.up_genes
    output:
        bar_chart = resultdir+"/GO/GO_{comp}_up.svg"
    message:
        "GO analysis of upregulated genes in {wildcards.comp} represented as a bar chart."
    conda:
        "gprofiler-1.0.0"
    resources:
        mem_mb = 3000,
        runtime = 60
    script:
        "../scripts/GO_bar_charts.py"


rule GO_downregulated_genes:
    input:
        genes = rules.volcano_plot.output.down_genes
    output:
        bar_chart = resultdir+"/GO/GO_{comp}_down.svg"
    message:
        "GO analysis of downregulated genes in {wildcards.comp} represented as a bar chart."
    conda:
        "gprofiler-1.0.0"
    resources:
        mem_mb = 3000,
        runtime = 60
    script:
        "../scripts/GO_bar_charts.py"
