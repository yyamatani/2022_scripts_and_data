## --- Filtering out genes with low variance --- ##

## R v3.6.0 (x64 bit) was used for running this code


## --- Step 1: Prepare data files -- ##
# Input file: Batch effect-corrected gene expression data among stem cell samples
args1 = commandArgs(trailingOnly=TRUE)[1]
input_file <- args1

# Output file: Batch effect-corrected gene expression data filtering low variance data
args2 = commandArgs(trailingOnly=TRUE)[2]
output_file <- args2


#=======================================================================#


## -- Step 2: Read data files -- ##
# Read input file:
input <- read.table(input_file, header=TRUE, row.names=1, sep="\t", quote="")
data <- as.matrix(input)


#=======================================================================#


## -- Step 3: Calculate coeffient of variation (CV) for removing low variance gene data
# 3-1. Calculate standard deviation (SD)
SD <- apply(data, 1, sd)
# 3-2. Calculate mean expression value of genes
mean_expression <- apply(data, 1, mean)
# 3-3. Calculate coeffient of variation (CV) for each gene
CV <- SD/mean_expression
data_with_CV <- cbind(data, CV)
# 3-4. Detect dataset with CV >= 0.8
data_with_higher_CV <- data_with_CV[CV >= 0.8, ]
filtered_data <- input[(rownames(data) %in% rownames(data_with_higher_CV)),]
# 3-5. Get data filtering low variance data as output data file
Genes <- rownames(filtered_data)
output <- cbind(Genes, filtered_data)
# 3-6. Get output data file
write.table(output, output_file, sep="\t", append=F, quote=F, row.names=F)


