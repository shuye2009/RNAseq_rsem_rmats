### Adapted from Danny Bergeron's code
resultdir = config['path']['resultdir']

rule gene_abundance_kallisto:
    input:
        quant = expand(rules.kallisto_quant.output.h5, sample=SAMPLES),
        samples = "config/design.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        out_file = resultdir+"/deseq2/kallisto_gene_level_abundance.tsv"
    params:
        kallisto_dir = resultdir+"/kallisto",
        out_dir = resultdir+"/deseq2",
        data_type = "kallisto"
    message:
        "Perform gene level abundance calculation for kallisto."
    conda:
        "deseq2-1.42.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/gene_level_abundance.R"

rule gene_abundance_rsem:
    input:
        quant = expand(rules.rsem_count.output.gene, sample=SAMPLES),
        samples = "config/design.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        out_file = resultdir+"/deseq2/rsem_gene_level_abundance.tsv"
    params:
        data_dir = resultdir+"/RSEM",
        out_dir = resultdir+"/deseq2",
        data_type = "rsem"
    message:
        "Perform gene level abundance calculation for RSEM."
    conda:
        "deseq2-1.42.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/gene_level_abundance.R"

rule deseq2_kallisto:
    input:
        quant = expand(rules.kallisto_quant.output.tsv, sample=SAMPLES),
        samples = "config/design.tsv",
        comparisons = "config/comparisons.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        out_files = resultdir+"/deseq2/{comp}_kallisto.csv",
    params:
        kallisto_dir = resultdir+"/kallisto",
        out_dir = resultdir+"/deseq2",
        data_type = "kallisto"
    message:
        "Perform differential expression analysis for {wildcards.comp}."
    conda:
        "deseq2-1.42.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/DESeq2_tximport.R"

rule deseq2_rsem:
    input:
        quant = expand(rules.rsem_count.output.isoform, sample=SAMPLES),
        samples = "config/design.tsv",
        comparisons = "config/comparisons.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        out_files = resultdir+"/deseq2/{comp}_rsem.csv",
    params:
        data_dir = resultdir+"/RSEM",
        out_dir = resultdir+"/deseq2",
        data_type = "rsem"
    message:
        "Perform differential expression analysis for {wildcards.comp} using RSEM quantification."
    conda:
        "deseq2-1.42.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/DESeq2_tximport.R"


rule volcano_plot_kallisto:
    input:
        DE_output = rules.deseq2_kallisto.output.out_files,
        filtered_genes = rules.merge_kallisto_quant.output.tpm,
    output:
        volcano = resultdir+"/deseq2/{comp}_kallisto.svg",
        up_genes = resultdir+"/deseq2/{comp}_kallisto_sig_DE_up.tsv",
        down_genes = resultdir+"/deseq2/{comp}_kallisto_sig_DE_down.tsv",
        all_genes = resultdir+"/deseq2/{comp}_kallisto_DE.tsv"
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

rule volcano_plot_rsem:
    input:
        DE_output = rules.deseq2_rsem.output.out_files,
        filtered_genes = rules.merge_rsem.output.tpm,
    output:
        volcano = resultdir+"/deseq2/{comp}_rsem.svg",
        up_genes = resultdir+"/deseq2/{comp}_rsem_sig_DE_up.tsv",
        down_genes = resultdir+"/deseq2/{comp}_rsem_sig_DE_down.tsv",
        all_genes = resultdir+"/deseq2/{comp}_rsem_DE.tsv"
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
        genes = rules.volcano_plot_rsem.output.up_genes
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
        genes = rules.volcano_plot_rsem.output.down_genes
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
