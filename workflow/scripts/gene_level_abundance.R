# https://bioconductor.org/packages/3.7/bioc/vignettes/tximport/inst/doc/tximport.html#use-with-downstream-bioconductor-dge-packages
### Script written by Danny Bergeron

#############################################################
#----------------- Loading the librairies -------------------
#############################################################
library(readr)
library(tximport)
library(DESeq2)

# Variables coming from the snakemake
data_dir <- snakemake@params[["data_dir"]]
output_dir <- snakemake@params[["out_dir"]]
design_file <- snakemake@input[["samples"]]
transcript_gene_file <- snakemake@input[["gene_id"]]
files <- snakemake@input[["quant"]]
gene_file <- snakemake@input[["gene_name"]]
data_type <- snakemake@params[["data_type"]]

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
sampleTable$condition <- factor(sampleTable$condition, levels=conditions)

# Read the transcript-gene matrix
tx2gene <- read_tsv(transcript_gene_file, col_names=c('TXNAME', 'GENEID'))

# Read the gene id-name matrix
id2name <- as.data.frame(read_tsv(gene_file, col_names=c('GENEID', 'GENENAME')))


# Get the h5 files for all conditions
#if(data_type == "kallisto"){
#    files <- file.path(data_dir, samples, "abundance.h5")
#}else if(data_type == "rsem"){
#    files <- file.path(data_dir, samples, paste0(samples,".isoforms.results"))
#}

names(files) <- samples
message("Gene level calculation")
txi.data <- tximport(files, type=data_type, tx2gene=tx2gene) # for gene differential analysis

# Output gene level abundance
message("Gene level abundance calculation")
gene_abundance <- as.data.frame(txi.data$abundance)
gene_abundance$GENEID <- row.names(gene_abundance)
gene_abundance <- merge(id2name, gene_abundance)
write_tsv(gene_abundance, file.path(output_dir, paste0(data_type, "_gene_level_abundance.tsv")))

# Output gene level counts
message("Gene level counts calculation")
gene_counts <- as.data.frame(txi.data$counts)
gene_counts$GENEID <- row.names(gene_counts)
gene_counts <- merge(id2name, gene_counts)
write_tsv(gene_counts, file.path(output_dir, paste0(data_type, "_gene_level_counts.tsv")))


# TO DO, use it to replace merge_kallisto and merge rsem
message("Isoform level calculation")
txi.data <- tximport(files, type=data_type, txIn=TRUE, txOut=TRUE) # for transcript differential analysis

# Output transcript level abundance 
message("Isoform level abundance calculation")
transcript_abundance <- as.data.frame(txi.data$abundance)
transcript_abundance$TXNAME <- row.names(transcript_abundance)
transcript_abundance <- merge(tx2gene, transcript_abundance)
write_tsv(transcript_abundance, file.path(output_dir, paste0(data_type, "_isoform_level_abundance.tsv")))

# Output transcript level counts
message("Isoform level counts calculation")
transcript_counts <- as.data.frame(txi.data$counts)
transcript_counts$TXNAME <- row.names(transcript_counts)
transcript_counts <- merge(tx2gene, transcript_counts)
write_tsv(transcript_counts, file.path(output_dir, paste0(data_type, "_isoform_level_counts.tsv")))

message("Gene and isoform level calculation completed")
