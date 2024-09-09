#!/bin/bash

source /cluster/home/t128737uhn/miniconda3/etc/profile.d/conda.sh
conda activate snakemake-7

snakemake \
    --unlock \
    --snakefile "/cluster/home/t128737uhn/snakemake_pipelines/RNAseq/workflow/Snakefile" \

snakemake \
    --cluster "sbatch -t {cluster.time} --mem {cluster.memory} -J {cluster.job-name} -p {cluster.partition} -N {cluster.nodes} -c {cluster.ntasks-per-node} -D {cluster.chdir} -o {cluster.output} -e {cluster.error}" \
    --cluster-config "/cluster/home/t128737uhn/snakemake_pipelines/RNAseq/profile_slurm/profile_workflow/cluster.json" \
    --snakefile "/cluster/home/t128737uhn/snakemake_pipelines/RNAseq/workflow/Snakefile" \
    --stats RNAseq_statistics.json \
    --latency-wait 2400 \
    --rerun-incomplete \
    --use-conda \
    --jobs 80


# 
conda deactivate