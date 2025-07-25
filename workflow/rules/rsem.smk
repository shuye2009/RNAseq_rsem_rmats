sampledir = config['path']['datadir']
resultdir = config['path']['resultdir']
indexdir = config["path"]["rsemdir"]
RSEM_SUFFIX = ["grp", "ti","chrlist","transcripts.fa", "seq","idx.fa","n2g.idx.fa"]
index_prefix = os.path.splitext(os.path.basename(config["path"]["genome_fasta"]))[0]

rule rsem_index:
    input:
        fasta = os.path.abspath(config["path"]["genome_fasta"]),
        gtf = os.path.abspath(config["path"]["genome_gtf"])
    output:
        rsem = expand(indexdir + "/" + index_prefix + ".{suffix}", suffix = RSEM_SUFFIX)
    params:
        rsem_prefix = index_prefix,
        rsem_dir = indexdir
    conda:
        "rsem-1.3.3"
    threads: 10
    resources:
        mem_mb = 10000,
        runtime = 240
    shell:
        """
        ln -f -s {input} {params.rsem_dir}/.
        rsem-prepare-reference -p {threads} --gtf {input.gtf} {input.fasta} {params.rsem_dir}/{params.rsem_prefix}
        """


rule rsem_count:
    input:
        index = expand(rules.rsem_index.output.rsem, suffix = RSEM_SUFFIX),
        bam = rules.star_align.output.txbam,
    output:
        gene = resultdir+"/RSEM/{sample}/{sample}.genes.results",
        isoform = resultdir+"/RSEM/{sample}/{sample}.isoforms.results"
    params:
        stranded = config["stranded"],
        prefixOut = resultdir+"/RSEM/{sample}/{sample}",
        rsem_index = indexdir + "/" + index_prefix
    log:
        resultdir+"/logs/RSEM/{sample}_rsem.log"
    message:
        "Quantitate transcript and gnee abundance for {wildcards.sample}."
    conda:
        "rsem-1.3.3"
    threads: 10
    resources:
        mem_mb = 10000,
        runtime = 240
    shell:
        "rsem-calculate-expression "
        "--no-bam-output "
        "--alignments "
        "--paired-end "
        "-p {threads} "
        "--strandedness {params.stranded} "
        "{input.bam} {params.rsem_index} {params.prefixOut}"

