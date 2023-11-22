args <- commandArgs(trailingOnly = TRUE)


orgdb <- args[1]
out_se <- args[2]
out_sce <- args[3]
#out_sizeFactors <- args[4]

starsolo_dir <- "results/STARsolo"
cogent_dir <- "results/cogent_analyze"
samplesheet <- "config/units.tsv"
kallisto_dir <- "results/kallisto_quant_tcc/"

# Packages loaded

library(dplyr)
library(stringr)
library(readr)
library(tibble)
# load the org.db for your organism
if(!require(orgdb, character.only=TRUE)){
    BiocManager::install(orgdb)
    library(orgdb, character.only=TRUE)
}
library(DESeq2)
library(AnnotationDbi)
library(tximport)
library(GenomicFeatures)
library(scater)
library(DropletUtils)

# Read STARsolo counts
files <- file.path(list.dirs(starsolo_dir, recursive=FALSE), "GeneFull/raw/")
names(files) <- str_remove_all(files, "\\S+STARsolo\\/|\\.Solo\\.out\\/GeneFull\\S+")
starsolo <- read10xCounts(files, names(files))
colnames(starsolo) <- starsolo$Sample

# Read TPMs

files <- file.path(list.dirs(kallisto_dir, recursive=FALSE), "total", "abundance.gene_1.tsv")
names(files) <- basename(str_remove(files, "\\/[^\\/]+\\/abundance.gene_1.tsv"))

kallisto_abund <- lapply(files, function(x){
    read_tsv(x, col_types="ccdd", col_names=TRUE)
})
kallisto_tpms <- Reduce(function(a,b) full_join(a, b, by="gene_id"), lapply(setNames(nm=names(kallisto_abund)), function(x){
    kallisto_abund[[x]] %>%
        dplyr::select(gene_id, tpm) %>%
        dplyr::rename({{x}} := "tpm")
})) %>% tibble::column_to_rownames("gene_id") %>% as.matrix()

kallisto_cts <- Reduce(function(a,b) full_join(a, b, by="gene_id"), lapply(setNames(nm=names(kallisto_abund)), function(x){
    kallisto_abund[[x]] %>%
        dplyr::select(gene_id, est_counts) %>%
        dplyr::rename({{x}} := "est_counts")
})) %>% tibble::column_to_rownames("gene_id") %>% as.matrix()

# Read CogentAP
files <- unlist(lapply(list.dirs(cogent_dir, recursive=FALSE), list.files, full.names=TRUE, pattern="_umi_uss_genematrix.csv"))
names(files) <- basename(str_remove(files, "_umi_uss_genematrix.csv"))

cogent_list <- lapply(setNames(nm=names(files)), function(x){
                    read_csv(files[[x]], col_names=c("gene_id", x), skip=1, col_types="cd")
                })
cogent <- Reduce(function(a, b) full_join(a, b, by="gene_id"), cogent_list)
cogent <- cogent %>% tibble::column_to_rownames("gene_id") %>% as.matrix()
rownames(cogent) <- str_extract(rownames(cogent), "^[^_]+")

# Row annot

# add gene symbols
gene_names_df <- data.frame(row.names = rownames(starsolo))

ens_no_version <- str_remove(rownames(gene_names_df), "\\.\\d+$")
stopifnot(length(ens_no_version) == length(unique(ens_no_version)))

gene_names_df$Symbol <- AnnotationDbi::mapIds(eval(as.name(orgdb)), ens_no_version, 
                                              keytype="ENSEMBL", column="SYMBOL", 
                                              multiVals="first")

gene_names_df$Uniq_syms <- scater::uniquifyFeatureNames(rownames(gene_names_df), gene_names_df$Symbol)
gene_names_df$entrez <- AnnotationDbi::mapIds(eval(as.name(orgdb)), ens_no_version, 
                                              keytype="ENSEMBL", column="ENTREZID", 
                                              multiVals="first") # there are duplicates in here.

gene_names_df$Gene_name <- AnnotationDbi::mapIds(eval(as.name(orgdb)), ens_no_version, 
                                                 keytype="ENSEMBL", column="GENENAME", 
                                                 multiVals="first")

# Sample meta

data_for_DE <- read_tsv(samplesheet) %>%
  as.data.frame() %>%
  dplyr::select(-fq1, -fq2) %>%
  unique() # samplesheet can have more than one row for a given sample (e.g. sequenced on more than one lane)

# samplesheet must have at least sample and group
stopifnot(c("sample", "group") %in% colnames(data_for_DE))

rownames(data_for_DE) <- data_for_DE$sample

# make sure order of samples in the meta data matches the counts
data_for_DE <- data_for_DE[colnames(starsolo), ]

stopifnot(all(!is.na(data_for_DE$group)))

# Make SingleCellExperiment

stopifnot(identical(sort(rownames(starsolo)), sort(rownames(cogent))))
stopifnot(identical(sort(rownames(starsolo)), sort(rownames(kallisto_tpms))))
stopifnot(identical(sort(rownames(starsolo)), sort(rownames(kallisto_cts))))

stopifnot(identical(sort(colnames(starsolo)), sort(colnames(cogent))))
stopifnot(identical(sort(colnames(starsolo)), sort(colnames(kallisto_tpms))))
stopifnot(identical(sort(colnames(starsolo)), sort(colnames(kallisto_cts))))

sce <- starsolo

assay(sce, "cogent") <- cogent[rownames(sce), colnames(sce)]
assay(sce, "tpms") <- kallisto_tpms[rownames(sce), colnames(sce)]
assay(sce, "kallisto_cts") <- kallisto_cts[rownames(sce), colnames(sce)]

rowData(sce) <- cbind(rowData(sce), gene_names_df[rownames(sce), ])
colData(sce) <- cbind(colData(sce), data_for_DE[colnames(sce), ])

sce <- sce[, order(sce$group)] # order samples by group


# Extract ERCC counts

sce2 <- splitAltExps(sce, ifelse(str_detect(rownames(sce), "^ERCC\\-\\d+$"), "ERCC", "gene"), ref = 'gene')

write_rds(sce2, out_sce)

se <- SummarizedExperiment(assays=list(counts=as.matrix(assay(sce2, "counts"))), colData=colData(sce2), rowData=rowData(sce2))

# Add vst
dds <- DESeqDataSet(se, design = ~ group)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

assays(se)$vst <- assay(vsd)

write_rds(se, out_se)



# Output size factors for use with other tools
#sizefactors <- 1/dds$sizeFactor
#tibble::tibble(sample=names(sizefactors), sizefactor=sizefactors) %>%
#  write_tsv(., out_sizeFactors) # "SizeFactors will now contain a factor for every sample which can be used to divide the 4th colun of a bedGraph/bigwig by. Both the aforementioned tools (bamCoverage and genomecov) have options though to directly scale the output by a factor (--scaleFactor or --scale respectively). !! Note though that these options will multiply rather than divide the 4th column with the factor, therefore you would need to provide the reciprocal as mentioned in the code chunk above." https://www.biostars.org/p/413626/#414440
