## --- Non-differentially expressed housekeeping genes (non-DEHKG) for RUVSeq --- ##

## R v3.6.0 (x64 bit) was used for running this code


## --- Step 1: Prepare data files -- ##
# Input file 1: Mean expression data of housekeeping genes
args1 = commandArgs(trailingOnly=TRUE)[1]
input_file_1 <- args1

# Input file 2: Expression data of genes that correspond with housekeeping genes
args2 = commandArgs(trailingOnly=TRUE)[2]
input_file_2 <- args2

# Output file: non-DEHKG expression data for RUVSeq
args3 = commandArgs(trailingOnly=TRUE)[3]
output_file <- args3


#=======================================================================#


## -- Step 2: Detect non-DEHKG data -- ##
# 2-1. Read input file 1: Mean expression data of housekeeping genes
input_1 <- read.table(input_file_1, header=TRUE, row.names=1, sep="\t", quote="")
data_1 <- as.matrix(input_1)

# 2-2. Calculate relative expression data and detect non-DEHKG data
mean_expssion_data <- apply(data_1, 1, mean)
calculate_relative_fold_change <- function(val_col) {return (log2((val_col/mean_expssion_data)+1))}
relative_log2_expression_data <- apply(data_1, 2, calculate_relative_fold_change)

# 2-3. Detect dataset with 0.5849625 < log2((val_col/mean_expssion_data)+1) < 1.584963 as non-DEHKG data
check_1 <- apply(relative_log2_expression_data, 1, function(x) all(x > 0.5849625))
non_DEHKG_list_candidate <- relative_log2_expression_data[check_1,]
check_2 <- apply(non_DEHKG_list_candidate, 1, function(x)  all(x < 1.584963))
non_DEHKG_list <- non_DEHKG_list_candidate[check_2,]
non_DEHKG_list <- as.matrix(non_DEHKG_list)


#=======================================================================#


## -- Step 3: Extract non-DEHKG expression data of 103 stem cell sample replicates -- ##
# 3-1. Read input file 2: Expression data of genes that correspond with housekeeping genes
input_2 <- read.table(input_file_2, header=TRUE, row.names=1, sep="\t", quote="")
data_2 <- as.matrix(input_2)

# 3-2. Extract non-DEHKG expression data of 103 stem cell sample replicates
non_DEHKG_expression_data <- data_2[(rownames(data_2) %in% rownames(non_DEHKG_list)),]
Genes <- rownames(non_DEHKG_expression_data)
output <- cbind(Genes, non_DEHKG_expression_data)
# 3-3. Get output data file
write.table(output, output_file, sep="\t", append=F, quote=F, row.names=F)


