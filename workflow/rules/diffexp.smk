### Adapted from Danny Bergeron's code
resultdir = config['path']['resultdir']

rule gene_abundance_kallisto:
    input:
        quant = expand(rules.kallisto_quant.output.h5, sample=SAMPLES),
        samples = "config/design.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        tpm_file = resultdir+"/deseq2/kallisto_gene_level_abundance.tsv",
        counts_file = resultdir+"/deseq2/kallisto_gene_level_counts.tsv",
        isoform_tpm_file = resultdir+"/deseq2/kallisto_isoform_level_abundance.tsv",
        isoform_counts_file = resultdir+"/deseq2/kallisto_isoform_level_counts.tsv"
    params:
        kallisto_dir = resultdir+"/kallisto",
        out_dir = resultdir+"/deseq2",
        data_type = "kallisto"
    message:
        "Perform gene and isoform level abundance and counts calculation for kallisto."
    conda:
        "deseq2-1.42.0"
    resources:
        mem_mb = 8000,
        runtime = 60
    script:
        "../scripts/gene_level_abundance.R"

rule gene_abundance_rsem:
    input:
        quant = expand(rules.rsem_count.output.isoform, sample=SAMPLES),
        samples = "config/design.tsv",
        gene_id = rules.tx2gene.output.tsv,
        gene_name = rules.tx2gene.output.gene_tsv
    output:
        tpm_file = resultdir+"/deseq2/rsem_gene_level_abundance.tsv",
        counts_file = resultdir+"/deseq2/rsem_gene_level_counts.tsv",
        isoform_tpm_file = resultdir+"/deseq2/rsem_isoform_level_abundance.tsv",
        isoform_counts_file = resultdir+"/deseq2/rsem_isoform_level_counts.tsv"
    params:
        data_dir = resultdir+"/RSEM",
        out_dir = resultdir+"/deseq2",
        data_type = "rsem"
    message:
        "Perform gene and isoform level abundance and counts calculation for RSEM."
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


