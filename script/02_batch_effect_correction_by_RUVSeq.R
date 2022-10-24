## --- Batch effect correction using RUVg function in RUVSeq (Risso et al., 2014) --- ##
## --- Using package RUVSeq version 1.18.0 at following Step 4 --- ##
## --- Remove Unwanted Variation from RNA-Seq Data (RUVSeq) is the remove unwanted variation (RUV) methods developed by Risso et al. (2014). -- ##
## --- For more information, please see: https://bioconductor.org/packages/release/bioc/html/RUVSeq.html --- ##

## R v3.6.0 (x64 bit) was used for running this code

## --- Step 1: Prepare data files -- ##
# Input file: Gene expression data among stem cell samples
args1 = commandArgs(trailingOnly=TRUE)[1]
input_file <- args1

# Reference file 1: House keeping gene data
args2 = commandArgs(trailingOnly=TRUE)[2]
reference_file <- args2

# Output file: Batch effect-corrected gene expression data among stem cell samples
args3 = commandArgs(trailingOnly=TRUE)[3]
output_file <- args3

# Reference file 2: Sample annotation data (e.g., stem cell type)
args4 = commandArgs(trailingOnly=TRUE)[4]
annotation_file <- args4


#=======================================================================#


## -- Step 2: Install and load package of RUVSeq version 1.18.0 (Risso et al., 2014) -- ##
# Install "RUVSeq" package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("RUVSeq", quietly = TRUE))
BiocManager::install("RUVSeq")

# Load "RUVSeq" package
library(RUVSeq)


#=======================================================================#


## -- Step 3: Read data files -- ##
# Read input file: Gene expression data among samples
input <- read.table(input_file, header=TRUE, row.names=1, sep="\t", quote="")
data <- as.matrix(input)

# Read reference file: House keeping gene data
reference <- read.table(reference_file, header=TRUE, row.names=1, sep="\t", quote="")
# Get House keeping genes 
house_keeping_gene <- rownames(reference)

# Read reference file: Sample annotation data (e.g., stem cell type)
annotation <- read.table(annotation_file, header=TRUE, sep="\t", quote="")
meta_data <- annotation
stem_cell_type <- c(meta_data$Stem_cell_type)


#=======================================================================#


## -- Step 4: Run RUVg function in RUVSeq version 1.18.0 (Risso et al., 2014) to estimate and correct batch effects -- ##
# 4-1. Get expression data of house keeping genes
cIdx <- rownames(data)[(rownames(data) %in% rownames(reference))]
# 4-2. Store data
set <- newSeqExpressionSet(data, phenoData = data.frame(stem_cell_type, row.names=colnames(data)))
# 4-3. Run RUVg: Estimating the factors of unwanted variation using control genes
set1 <- RUVg(set, cIdx, k=1)
# 4-4. Get Expression data table after batch effect correction
batch_effect_corrected_data <- normCounts(set1)
# 4-5. Remove read count data with less than three counts in all stem cell samples
batch_effect_corrected_data[batch_effect_corrected_data < 3] <- 0
batch_effect_corrected_data <- batch_effect_corrected_data[apply(batch_effect_corrected_data,1,sum)>0,]


#=======================================================================#


## -- Step 5: Output data file after batch effect correction -- ##
# 5-1. Get batch effect-corrected expression data as output data file
Genes <- rownames(batch_effect_corrected_data)
output <- cbind(Genes, batch_effect_corrected_data)
# 5-2. Get output data file
write.table(output, file=output_file, sep = "\t", quote=FALSE, row.names=FALSE)


