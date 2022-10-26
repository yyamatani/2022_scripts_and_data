## --- Calculating the ratio of between-group and within-group variance before and after batch effect correction in research groups --- ##

## R v3.6.0 (x64 bit) was used for running this code


## --- Step 1: Prepare data files -- ##
# Input file 1: raw gene expression data or batch effect-corrected gene expression data
args1 = commandArgs(trailingOnly=TRUE)[1]
input_file_raw <- args1

# Input file 2: Batch effect-corrected gene expression data among stem cell samples
args2 = commandArgs(trailingOnly=TRUE)[2]
input_file_corrected <- args2


#=======================================================================#


## -- Step 2: Calculating the ratio of between-group and within-group variance before and after the correction in research groups -- ##
before_or_after_list <- c("before", "after")
data_list <- c(input_file_raw, input_file_corrected)

for (i in 1:2) {
# 3-1. Read input file:
current_data <- read.table(data_list[i], header=TRUE, row.names=1, sep="\t", quote="")
data <- as.matrix(current_data)

# 3-2. Extract dataset of each research group (Nos.1-24)
No1_data <- data[,c(1:2)]
No2_data <- data[,c(3:4)]
No3_data <- data[,c(5:7)]
No4_data <- data[,c(8:9)]
No5_data <- data[,c(10:12)]
No6_data <- data[,c(13:14)]
No7_data <- data[,c(15:16)]
No8_data <- data[,c(17:19)]
No9_data <- data[,c(20:26)]
No10_data <- data[,c(27:28)]
No11_data <- data[,c(29:33)]
No12_data <- data[,c(34:35)]
No13_data <- data[,c(36:38)]
No14_data <- data[,c(39:41)]
No15_data <- data[,c(42:59)]
No16_data <- data[,c(60:63)]
No17_data <- data[,c(64:65)]
No18_data <- data[,c(66:69)]
No19_data <- data[,c(70:80)]
No20_data <- data[,c(81:84)]
No21_data <- data[,c(85:93)]
No22_data <- data[,c(94:96)]
No23_data <- data[,c(97:98)]
No24_data <- data[,c(99:103)]

# 3-3. Calculate mean value of research group data replicates
number_of_research_groups <- 24

for (j in 1:number_of_research_groups) {
No_data <- eval(parse(text = paste("No", j, "_data", sep = "")))
assign(paste("No", j, "_mean_expression", sep=""), colMeans(No_data))
}
research_group_expression_data <- cbind(No1_mean_expression, No2_mean_expression, No3_mean_expression, No4_mean_expression,
					No5_mean_expression, No6_mean_expression, No7_mean_expression, No8_mean_expression, 
					No9_mean_expression, No10_mean_expression, No11_mean_expression, No12_mean_expression, No13_mean_expression, No14_mean_expression, No15_mean_expression, No16_mean_expression, No17_mean_expression, No18_mean_expression, No19_mean_expression,
					No20_mean_expression, No21_mean_expression, No22_mean_expression, No23_mean_expression, No24_mean_expression
					)

# 3-4. Calculate within-group variance
sample_mean_expression <- colMeans(research_group_expression_data)
sample_mean_expression_matrix <- matrix(rep(sample_mean_expression, nrow(research_group_expression_data)), nrow=nrow(research_group_expression_data), ncol=ncol(research_group_expression_data),byrow=TRUE)
within_group <- research_group_expression_data - sample_mean_expression_matrix
within_group_squared <- sum(within_group^2)
within_group_variance <- within_group_squared/((nrow(research_group_expression_data)-1)*ncol(research_group_expression_data))

# 3-5. Calculate between-group variance
overall_mean_expression <- sum(colMeans(research_group_expression_data))/ncol(research_group_expression_data)
between_group <- matrix(rep(sample_mean_expression - overall_mean_expression, nrow(research_group_expression_data)),nrow(research_group_expression_data), byrow=TRUE)
between_group_squared <- sum(between_group^2)
between_group_variance <- between_group_squared/(ncol(research_group_expression_data)-1)

# 3-6. Calculate the ratio of variance and print the result 
X <- paste("###### The ratio of variance ", before_or_after_list[i], " batch effect correction in research groups ###### ", sep="")
print(X)
print("between-group variance in research groups")
print(between_group_variance)
print("within-group variance in research groups")
print(within_group_variance)
print("between-group variance / within-group variance = the ratio of variance in research group")
print(between_group_variance/within_group_variance)
}

