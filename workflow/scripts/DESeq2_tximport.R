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
comparison_file <- snakemake@input[["comparisons"]]
transcript_gene_file <- snakemake@input[["gene_id"]]
gene_file <- snakemake@input[["gene_name"]]
files <- snakemake@input[["quant"]]
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

# Loading the comparisons
comparisons <- read.table(
    comparison_file, header=TRUE
)

# Read the transcript-gene matrix
tx2gene <- read_tsv(transcript_gene_file, col_names=c('TXNAME', 'GENEID'))

# Read the gene id-name matrix
id2name <- as.data.frame(read_tsv(gene_file, col_names=c('GENEID', 'GENENAME')))


# Get the h5 files for all conditions
#if(data_type == "kalliston"){
#    files <- file.path(data_dir, samples, "abundance.h5")
#}else if(data_type == "rsem"){
#    files <- file.path(data_dir, samples, paste0(samples,".isoforms.results"))
#}

names(files) <- samples
txi.data <- tximport(files, type=data_type, tx2gene=tx2gene) # for gene differential analysis
# txi.data <- tximport(files, type=data_type, txOut=TRUE) # for transcript differential analysis


#############################################################
#--------------- Creating the DESeq2 object -----------------
#############################################################
# Create the DESeq2 object with all the conditions
dds <- DESeqDataSetFromTximport(
    txi.data,
    sampleTable,
    ~condition
)

# pre filtering the count
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# put the levels for each columns
dds$condition <- factor(dds$condition, conditions)

# call the function for differential gene expression analysis
dds <- DESeq(dds)


#############################################################
#-------------- Looping over all comparisons ----------------
#############################################################
for (row in 1:nrow(comparisons)) {

    condition1 <- comparisons[row, "cdn1"]
    condition2  <- comparisons[row, "cdn2"]
    exp <- sprintf("%s-%s", condition1, condition2)

    res <- results(
        dds,
        contrast = c(
            "condition",
            condition2,
            condition1
        ),
        #pAdjustMethod = "fdr"
    )

    # transform the result to data frame
    res_df <- as.data.frame(res)

    # Writing results to file
    fname <- paste(
        output_dir,
        paste(exp, "csv", sep='.'),
        sep='/'
    )

    write.csv(
        res_df,
        file=fname,
        quote=FALSE,
    )
}