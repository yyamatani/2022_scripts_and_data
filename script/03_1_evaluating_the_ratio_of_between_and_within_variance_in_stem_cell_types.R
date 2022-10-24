## --- Calculating the ratio of between-group and within-group variance before and after batch effect correction in stem cell types --- ##

## R v3.6.0 (x64 bit) was used for running this code


## --- Step 1: Prepare data files -- ##
# Input file 1: raw gene expression data or batch effect-corrected gene expression data
args1 = commandArgs(trailingOnly=TRUE)[1]
input_file_raw <- args1

# Input file 2: Batch effect-corrected gene expression data among stem cell samples
args2 = commandArgs(trailingOnly=TRUE)[2]
input_file_corrected <- args2


#=======================================================================#


## -- Step 2: Calculating the ratio of between-group and within-group variance before and after the correction in stem cell types -- ##
before_or_after_list <- c("before", "after")
data_list <- c(input_file_raw, input_file_corrected)

for (i in 1:2) {
# 3-1. Read input file:
current_data <- read.table(data_list[i], header=TRUE, row.names=1, sep="\t", quote="")
data <- as.matrix(current_data)

# 3-2. Extract dataset of each stem cell type
ESC_data <- data[,c(1:9)] # ESC
iPSC_data <- data[,c(10:19)] # iPSC
MSC_data <- data[,c(20:80)] # MSC
HSC_data <- data[,c(81:103)] # HSC

# 3-3. Calculate mean value of stem cell sample replicates
ESC_mean_expression <- colMeans(ESC_data)
iPSC_mean_expression <- colMeans(iPSC_data)
MSC_mean_expression <- colMeans(MSC_data)
HSC_mean_expression <- colMeans(HSC_data)
stem_cell_expression_data <- cbind(ESC_mean_expression, iPSC_mean_expression, MSC_mean_expression, HSC_mean_expression)

# 3-4. Calculate within-group variance
sample_mean_expression <- colMeans(stem_cell_expression_data)
sample_mean_expression_matrix <- matrix(rep(sample_mean_expression, nrow(stem_cell_expression_data)), nrow=nrow(stem_cell_expression_data), ncol=ncol(stem_cell_expression_data),byrow=TRUE)
within_group <- stem_cell_expression_data - sample_mean_expression_matrix
within_group_squared <- sum(within_group^2)
within_group_variance <- within_group_squared/((nrow(stem_cell_expression_data)-1)*ncol(stem_cell_expression_data))

# 3-5. Calculate between-group variance
overall_mean_expression <- sum(colMeans(stem_cell_expression_data))/ncol(stem_cell_expression_data)
between_group <- matrix(rep(sample_mean_expression - overall_mean_expression, nrow(stem_cell_expression_data)),nrow(stem_cell_expression_data), byrow=TRUE)
between_group_squared <- sum(between_group^2)
between_group_variance <- between_group_squared/(ncol(stem_cell_expression_data)-1)

# 3-6. Calculate the ratio of variance and print the result 
X <- paste("###### The ratio of variance ", before_or_after_list[i], " batch effect correction in stem cell types ###### ", sep="")
print(X)
print("between-group variance in stem cell types")
print(between_group_variance)
print("within-group variance in stem cell types")
print(within_group_variance)
print("between-group variance / within-group variance = the ratio of variance in stem cell types")
print(between_group_variance/within_group_variance)
}

