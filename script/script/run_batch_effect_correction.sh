#!/bin/bash
#$ -S /bin/bash


# Set environment variables in R 3.6
module load R/3.6

root_dir=`pwd`"/" # Get root directory
script_dir=${root_dir}"script/" # Script directory
input_dir=${root_dir}"input/" # Input directory
reference_dir=${root_dir}"reference/" # Reference dir
output_dir=${root_dir}"output/"; if [ ! -e ${output_dir} ]; then mkdir -p ${output_dir}; fi # If output directory is not exist, make it

# Prepare scripts
script_01_nonDEHKG_detection_for_RUVSeq=${script_dir}"01_nonDEHKG_detection_for_RUVSeq.R"
script_02_batch_effect_correction_by_RUVSeq=${script_dir}"02_batch_effect_correction_by_RUVSeq.R"
script_03_1_evaluating_the_ratio_of_between_and_within_variance_in_stem_cell_types=${script_dir}"03_1_evaluating_the_ratio_of_between_and_within_variance_in_stem_cell_types.R"
script_03_2_evaluating_the_ratio_of_between_and_within_variance_in_research_groups=${script_dir}"03_2_evaluating_the_ratio_of_between_and_within_variance_in_research_groups.R"
script_04_filtering_out_genes_with_low_variance=${script_dir}"04_filtering_out_genes_with_low_variance.R"


#=================#


# Prepare raw expression data table for the evaluation of batch effect correction
raw_data_table=${input_dir}"raw_read_count_table_before_RUVSeq.txt"

# Prepare TPM normalized gene expression data table as input file of RUVSeq (TPM: Transcripts Per Million)
unzip ${input_dir}"tpm_read_count_table_before_RUVSeq.txt.zip" # Get "tpm_read_count_table_before_RUVSeq.txt"
input_data_table=${input_dir}"tpm_read_count_table_before_RUVSeq.txt"

# Prepare sample annotation list
sample_annotation_list=${reference_dir}"stem_cell_sample_annotation.txt"

# Prepare human housekeeping gene list downlowed from http://www.tau.ac.il/~elieis/HKG/
housekeeping_gene_list=${reference_dir}"housekeeping_gene_list.txt"

# Prepare batch-effect corrected gene expression data table as output file
output_data_table_after_RUVSeq=${output_dir}"read_count_table_after_RUVSeq.txt"
if [ ! -e ${output_data_table_after_RUVSeq} ]; then touch ${output_data_table_after_RUVSeq}; fi # If output file is not exist, make it

# Prepare batch-effect corrected gene expression data table with filtering low variance data as the final output file
output_data_table_after_RUVSeq_filtering_low_variance_data=${output_dir}"read_count_table_after_RUVSeq_filtering_low_variance_data.txt"
if [ ! -e ${output_data_table_after_RUVSeq_filtering_low_variance_data} ]; then touch ${output_data_table_after_RUVSeq_filtering_low_variance_data}; fi # If output file is not exist, make it


#=================#


## Step 1:
# For using housekeeping gene expression profile as negative control of batch effect correction by RUVSeq, select non-differentially expressed house keepig genes (non-DEHKG) following four steps

# 1-1. Make temporary directory 1
temporary_dir_1=${output_dir}"01_temporary_dir/"; if [ ! -e ${temporary_dir_1} ]; then mkdir -p ${temporary_dir_1}; fi


# 1-2. Calculating mean value of each stem cell sample
sample_replicate_list=(9 3 2 5 12 30 9 5 5 12 8 3) # List of replicates among 12 kinds of stem cell samples
number_of_replicates=${#sample_replicate_list[@]}; number_of_loop=$((number_of_replicates-1))

# Extract expression data of genes that correspond with housekeeping list genes
housekeeping_gene_data_table=${output_dir}"expression_data_table_of_genes_correspond_with_housekeeping_list_genes.txt" # Prepare expression data table of genes matched to housekeeping genes
temporary_file_1_1=${temporary_dir_1}"1_1_temporary_file.txt"
sed -e '1d' ${input_data_table} | cut -f1 | awk -F "\t" '{sub("_dup[0-9]",""); OFS="\t"; print $0;}' | sort | uniq > $temporary_file_1_1
cut -f1 ${housekeeping_gene_list} | sort | uniq | fgrep -xf "-" ${temporary_file_1_1} | while read line
do
awk -F "\t" -v check_housekeeping_gene=${line} '{sub("_dup[0-9]","",$1); OFS="\t"; if($1==check_housekeeping_gene){print $0;}}' ${input_data_table} | \
awk -F "\t" '{i=1; while(i <= NR) {gene=$1"_dup"i; i++;} OFS="\t"; print gene,$0;}' | cut -f1,3- >> ${housekeeping_gene_data_table}
done
unlink ${temporary_file_1_1}
header_row=`head -n1 ${input_data_table}`; sed -i "1s/^/${header_row}\n/" ${housekeeping_gene_data_table} # Add header row

for n in `seq 0 ${number_of_loop}`
do
number_of_samples=$((n+1))
if [ ${number_of_samples} -lt 10 ]; then number_of_samples="0"${number_of_samples}; fi
temporary_file_1_2=${temporary_dir_1}"1_2_temporary_file_stem_cell_sample"${number_of_samples}".txt"; if [ $n -eq 0 ]; then first_column=2; else first_column=$((end_column+1)); fi
end_column=$((first_column+sample_replicate_list[$n]-1))
sed -e '1d' ${housekeeping_gene_data_table} | cut -f${first_column}-${end_column} | \
awk -F "\t" '{{i=1;sum=0; while(i <= NF) {sum += $i; i++;} sample_mean=sum/NF; OFS="\t"; print sample_mean;}}' > ${temporary_file_1_2}
done
 
paste_file_list=`ls -1 ${temporary_dir_1} | sed -e 's/\s/\t/g' | tr "\n" " " | sed -e 's/\s$/\n/g' | less`
cd ${temporary_dir_1}
mean_expression_data_of_housekeeping_genes=${output_dir}"mean_expression_data_table_of_housekeeping_genes.txt"
sed -e '1d' ${housekeeping_gene_data_table} | cut -f1 | paste -d "\t" "-" ${paste_file_list} | \
sed "1s/^/Genes\tESC\tAMiPSC\tPBiPSC\tSFiPSC\tADMSC\tBMMSC\tUCMSC\tAMMSC\tCPMSC\tUCHSC\tBMHSC\tFLHSC\n/" > ${mean_expression_data_of_housekeeping_genes}
cd ${output_dir}


# 1-3. Remove temporary directory 1
rm -rf ${temporary_dir_1}


# 1-4. Detect non-DEHKG among 12 kinds of stem cell samples
if [ -z ${mean_expression_data_of_housekeeping_genes} ]; then echo "Please prepare mean expression data of housekeeping genes!"; exit ; fi # check file existence
if [ -z ${housekeeping_gene_data_table} ]; then echo "Please prepare expression data table of genes matched to housekeeping genes!"; exit ; fi # check file existence

housekeeping_genes_for_ruvseq=${output_dir}"housekeeping_gene_data_table_for_RUVSeq.txt" # prepare housekeeping gene expression data for RUVSeq
Rscript ${script_01_nonDEHKG_detection_for_RUVSeq} ${mean_expression_data_of_housekeeping_genes} ${housekeeping_gene_data_table} ${housekeeping_genes_for_ruvseq}


#=================#


## Step 2:
# Conduct batch effect correction by RUVSeq

# 2-1. Check file existence
if [ -z ${input_data_table} ]; then echo "Please prepare input data table!"; exit ; fi
if [ -z ${housekeeping_genes_for_ruvseq} ]; then echo "Please prepare housekeepig gene expression data for RUVSeq!"; exit ; fi
if [ -z ${output_data_table_after_RUVSeq} ]; then echo "Please prepare batch-effect corrected gene expression data table as output file!"; exit ; fi
if [ -z ${sample_annotation_list} ]; then echo "Please prepare sample annotation list!"; exit ; fi

# 2-2. Convert input data to integer and run batch effect correction by RUVSeq
awk -F "\t" '{if(NR==1){print $0;}}; NR>1{i=2; while(i <= NF) {$i=int($i); i++;} OFS="\t"; print $0;}' ${input_data_table} | \
Rscript ${script_02_batch_effect_correction_by_RUVSeq} stdin ${housekeeping_genes_for_ruvseq} ${output_data_table_after_RUVSeq} ${sample_annotation_list}
unlink ${housekeeping_gene_data_table}; unlink ${mean_expression_data_of_housekeeping_genes}; unlink ${housekeeping_genes_for_ruvseq}


#=================#


## Step 3: Calculate the ratio of between-group and within-group variance and compare it before and after batch effect correction in stem cell types and research groups, respectively
# 3-1. Calculate the ratio of between-group and within-group variance and print the result in stem cell types
Rscript ${script_03_1_evaluating_the_ratio_of_between_and_within_variance_in_stem_cell_types} ${raw_data_table} ${output_data_table_after_RUVSeq}
# 3-2. Calculate the ratio of between-group and within-group variance and print the result in research groups
Rscript ${script_03_2_evaluating_the_ratio_of_between_and_within_variance_in_research_groups} ${raw_data_table} ${output_data_table_after_RUVSeq}

# If the ratio of variance (between-group variance / within-group variance) become smaller after the correction than before, the within-group variance is larger after the correction.


#=================#


## Step 4: Filter out genes with low variance

# 4-1. Make temporary directory 2
temporary_dir_2=${output_dir}"02_temporary_dir/"; if [ ! -e ${temporary_dir_2} ]; then mkdir -p ${temporary_dir_2}; fi

# 4-2. Merge multiple isoform data by calculating mean value of them 
temporary_file_2_1=${temporary_dir_2}"2_1_temporary_file.txt"; if [ ! -e ${temporary_file_2_1} ]; then touch ${temporary_file_2_1}; fi
cut -f1 ${output_data_table_after_RUVSeq} | sed -e '1d' | awk -F "\t" '{sub("_dup[0-9]",""); print $0;}' | sort | uniq | while read line
do
sed -e '1d' ${output_data_table_after_RUVSeq} | awk -F "\t" '{sub("_dup[0-9]","",$1); OFS="\t"; print $0;}' | \
awk -F  "\t" -v current_unique_gene=${line} '{if($1==current_unique_gene){OFS="\t"; print $0;}}' | \
awk '{for(i=2;i<=NF;i++) total[i]+=$i;} END{for(i=2;i<=NF;i++) printf "\t%d ",total[i]/NR;}' | awk -F  "\t" '{OFS="\t"; print $0;}' >> ${temporary_file_2_1}
done
temporary_file_2_2=${temporary_dir_2}"2_2_temporary_file.txt"
cut -f1 ${output_data_table_after_RUVSeq} | sed -e '1d' | awk -F "\t" '{sub("_dup[0-9]",""); print $0;}' | sort | uniq | paste -d "" "-" ${temporary_file_2_1} > ${temporary_file_2_2}
temporary_file_2_3=${temporary_dir_2}"2_3_temporary_file.txt"
head -n1 ${output_data_table_after_RUVSeq} | cat "-" ${temporary_file_2_2} > ${temporary_file_2_3}

# 4-4. Filter out genes with low variance
Rscript ${script_04_filtering_out_genes_with_low_variance} ${temporary_file_2_3} ${output_data_table_after_RUVSeq_filtering_low_variance_data}

# 4-5. Remove temporary directory 2
rm -rf ${temporary_dir_2}


exit


