#!/sur/bin/Rscript
setwd("./HDAMS")

args <- commandArgs()
combine_normalized_count_file <-args[6]
combine_ratio_file <-args[7]
predicted_seq_specificity <-args[8]


normalized_counts_data <-read.table(combine_normalized_count_file, header =TRUE)

ma6_ratio_data <-read.table(combine_ratio_file, header =TRUE)

beta_level_data <- read.table(predicted_seq_specificity, header=TRUE)


combined_processed_data <-merge( normalized_counts_data, ma6_ratio_data,
                                 by.x="Id", by.y="Id", incomparables =NA )

combined_processed_data <-merge(combined_processed_data, beta_level_data, by.x="Id", by.y="seq_id", incomparable= NA)

combined_cols <-colnames(combined_processed_data)

ctr_ip_cols <-combined_cols[grep('Ip[0-9]*_ctr' ,combined_cols)]
ctr_input_cols <-combined_cols[grep('Input[0-9]*_ctr' ,combined_cols)]

trt_ip_cols <-combined_cols[grep('Ip[0-9]*_trt' ,combined_cols)]
trt_input_cols <-combined_cols[grep('Input[0-9]*_trt' ,combined_cols)]

combined_ip_cols <-append(ctr_ip_cols, trt_ip_cols)
combined_input_cols <-append(ctr_input_cols, trt_input_cols)

ctr_m6a_ratio <-combined_cols[grep('m6a_ratio_ctr[0-9]*' ,combined_cols)]
trt_m6a_ratio <-combined_cols[grep('m6a_ratio_trt[0-9]*' ,combined_cols)]
combined_m6a_ratio <-append(ctr_m6a_ratio, trt_m6a_ratio)

seq_beta_id <-combined_cols[grep('prob_beta_.*' ,combined_cols)]
  
S1_row <-dim(combined_processed_data)[1]

flag <-rep(TRUE, times=S1_row)

for (i in seq(1:S1_row))
{
  if(sum(combined_processed_data[i, combined_input_cols] <1) >0)
   
    { flag[i] <-FALSE }
}

combined_processed_data <-combined_processed_data[flag, ]
# cat(dim(combined_processed_data))

fout_connet <-file('Control_Treatment_Normalized_Count_diff.txt','w')

head_line <-paste('Site_id','Ip_counts','Input_counts','m6a_ratio','Sample_group_id','Beta_level', sep= "\t")

writeLines(head_line, fout_connet)

row_num <-dim(combined_processed_data)[1]

for (i in seq(1:row_num))
{ 
  ls_ip_count <-paste(combined_processed_data[i, combined_ip_cols], collapse=';')
  
  ls_input_count <-paste(combined_processed_data[i, combined_input_cols], collapse=';')
  
  ls_m6a_ratio <-paste(round(combined_processed_data[i, combined_m6a_ratio] ,3), collapse=';')
  
  ls_group_id <-paste(append(rep(1, length(ctr_m6a_ratio)), rep(2, length(trt_m6a_ratio))), collapse =';')
  
  ls_beta_id <-paste(round(combined_processed_data[i, seq_beta_id] ,3), collapse=';')
  
  collection_vec <-paste(combined_processed_data[i, 'Id'], ls_ip_count, 
                         ls_input_count, ls_m6a_ratio, ls_group_id, ls_beta_id, sep="\t")
  
  writeLines(collection_vec, fout_connet)
}

close(fout_connet)



