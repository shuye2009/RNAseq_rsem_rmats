# https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html#use-with-downstream-bioconductor-dge-packages
### Script written by Danny Bergeron
### Modified by Shuye Pu, used only when the number of samples in each condition is less than 2, thus deseq2 cannot be called

#############################################################
#----------------- Loading the librairies -------------------
#############################################################
library(readr)
library(tximport)
library(DESeq2)

# Variables coming from the snakemake
kallisto_dir <- snakemake@params[["kallisto_dir"]]
output_dir <- snakemake@params[["out_dir"]]
design_file <- snakemake@input[["samples"]]
transcript_gene_file <- snakemake@input[["gene_id"]]
gene_file <- snakemake@input[["gene_name"]]
out_file <- snakemake@output[["out_file"]]
# create the directory that will contain the results
dir.create(output_dir, showWarnings=FALSE)

#############################################################
#------------- Importing data and information ---------------
#############################################################
# Loading samples information from the design file.
sampleTable <- read.table(
    design_file, header=TRUE,
    row.names="sample", check.names=FALSE
)
conditions <- unique(sampleTable$condition)
samples <- rownames(sampleTable)

# Read the transcript-gene matrix
tx2gene <- read_tsv(transcript_gene_file, col_names=c('TXNAME', 'GENEID'))

# Read the gene id-name matrix
id2name <- as.data.frame(read_tsv(gene_file, col_names=c('GENEID', 'GENENAME')))

# Get the h5 files for all conditions
files <- file.path(kallisto_dir, samples, "abundance.h5")
names(files) <- samples
txi.kallisto <- tximport(files, type="kallisto", tx2gene=tx2gene) # for gene differential analysis

# Output gene level abundance
gene_abundance <- txi.kallisto$abundance
gene_id <- data.frame(GENEID = row.names(gene_abundance))
gene_name <- merge(id2name, gene_id)
gene_abundance <- cbind(gene_name, gene_abundance)
write_tsv(gene_abundance, out_file)
