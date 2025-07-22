resultdir = config['path']['resultdir']



rule multiqc:
    input:
        fastqc_pretrim = expand(rules.fastqc_pretrim.output,sample=SAMPLES,read=['R1','R2']),
        fastqc_posttrim = expand(rules.fastqc_posttrim.output.html,sample=SAMPLES,read=['R1','R2']),
        picard = expand(rules.picard.output.pdf,sample=SAMPLES),
        kallisto = expand(rules.kallisto_quant.output,sample=SAMPLES),
        rsem = expand(rules.rsem_count.output,sample=SAMPLES)
    output:
        html = resultdir+"/multiqc/multiqc_report.html"
    params:
        scan_dir = f"{resultdir}/fastqc/pretrim {resultdir}/fastqc/posttrim {resultdir}/STAR/*/*Log.final.out {resultdir}/logs/kallisto {resultdir}/logs/RSEM {resultdir}/picard {resultdir}/logs/trimmomatic",
        outdir = directory(resultdir+"/multiqc")
    log:
        resultdir+"/logs/multiqc.log"
    message:
        "Summarize analysis results for multiple tools and samples in a single report."
    conda:
        "fastqc-0.12.1"
    resources:
        mem_mb = 3000,
        runtime = 60
    shell:
        """
        if [[ -d {params.outdir}/multiqc_data ]]; then rm -r {params.outdir}/multiqc_data; fi
        multiqc {params.scan_dir} --outdir {params.outdir} &> {log}
        
        """
