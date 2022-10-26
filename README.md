Comprehensive Comparison of Gene Expression Diversity among a Variety of Human Stem Cells

Our study aimed to rank the relative importance of biological and technical factors to gene expression among undifferentiated human stem cell types using publicly available RNA-seq datasets.
To analyze RNA-seq data of undifferentiated human stem cells reported by different research groups,
we corrected batch effects in the gene expression data.
The batch effect correction was conducted using "Remove Unwanted Variation from RNA-Seq Data" (RUVSeq) developed by Risso et al. (2014).


sh script/run_batch_effect_correction.sh


The executable code, "run_batch_effect_correction.sh" runs batch effect correction following steps:

Step 1: For using the housekeeping gene expression profile as a negative control of batch effect correction by RUVSeq, select non-differentially expressed housekeeping genes (non-DEHKG) following four steps

Step 2: Conduct batch effect correction by RUVSeq version 1.18.0 (Risso et al., 2014)

Step 3: Calculate the ratio of between-group and within-group variance and compare it before and after batch effect correction in stem cell types and research groups, respectively

Step 4: Filter out genes with low variance
