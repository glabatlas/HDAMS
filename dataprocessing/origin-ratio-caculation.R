#!/sur/bin/Rscript
setwd("./HDAMS")

args <- commandArgs()
ctr_file <-args[6]
contr_group_num <-as.integer(args[7])

trt_file <-args[8]
treat_group_num <-as.integer(args[9])

contr_data <-read.table(ctr_file, header = TRUE)
treat_data <-read.table(trt_file, header = TRUE)


#########################################################################################################
###    Control group processing 
#########################################################################################################

event_id_ctr <- paste(contr_data[ ,"peakname"], contr_data[ ,"chrom"], contr_data[ ,"strand"], contr_data[ ,"start"], contr_data[ ,"end"], sep="_")

contr_data <- cbind(Id=event_id_ctr, contr_data)

S1 <- dim(contr_data)[1]

flag_1 <- rep(TRUE, times=S1)

ctr_col_name <- colnames(contr_data)

ctr_count_columns <- ctr_col_name[grep('Input[0-9]*$|Ip[0-9]*$',ctr_col_name)]

ctr_ip_count_columns <-ctr_col_name[grep('Ip[0-9]*$',ctr_col_name)]

for (i in seq(1:S1))
{   
  if (sum(c(contr_data[i, ctr_count_columns])!=0)< 2*contr_group_num | mean(as.numeric(contr_data[i, ctr_ip_count_columns])) <3)
    
  { flag_1[i] <- FALSE }
}

contr_data <- contr_data[flag_1, ]

combined_ctr_ratio_data <- contr_data

for (i in seq(1:contr_group_num))
{
  ip_id <- paste('Ip', as.character(i), sep="")
  ip_total_id <-paste('Ip', as.character(i), '_tatal_counts',sep="")
  input_id <-paste('Input', as.character(i), sep="")
  input_total_id <-paste('Input', as.character(i), '_tatal_counts', sep="")
  ctr_m6a_ratio <-paste('m6a_ratio_ctr_', as.character(i), sep="")

  ctr_ratio <- (as.numeric(contr_data[ , ip_id])/as.numeric(contr_data[ , ip_total_id]))/ (as.numeric(contr_data[ , input_id])/as.numeric(contr_data[ , input_total_id]))
  
  temp_col<- append(colnames(combined_ctr_ratio_data), ctr_m6a_ratio)
  combined_ctr_ratio_data <- cbind(combined_ctr_ratio_data, ctr_ratio)
  colnames(combined_ctr_ratio_data) <- temp_col
}

# write.table(combined_ctr_ratio_data, file="Control_combined_counts_ratio_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#########################################################################################################
###    Treatment group processing 
#########################################################################################################

event_id_trt <- paste(treat_data[ ,"peakname"], treat_data[ ,"chrom"], treat_data[ ,"strand"], treat_data[ ,"start"], treat_data[ ,"end"], sep="_")

treat_data <- cbind(Id=event_id_trt, treat_data)

S2 <- dim(treat_data)[1]

flag_2 <- rep(TRUE, times=S2)

trt_col_name <- colnames(treat_data)
trt_count_columns <- trt_col_name[grep('Input[0-9]*$|Ip[0-9]*$', trt_col_name)]
trt_ip_count_columns <-trt_col_name[grep('Ip[0-9]*$', trt_col_name)]

for (i in seq(1:S2))
{   
  if ( sum(c(treat_data[i, trt_count_columns])!=0)< 2*treat_group_num | mean(as.numeric(treat_data[i, trt_ip_count_columns]))<3 )
  
    { flag_2[i]<-FALSE }
}

treat_data <- treat_data[flag_2, ]

combined_trt_ratio_data <-treat_data

for (i in seq(1:treat_group_num))
{
  ip_id <- paste('Ip', as.character(i), sep="")
  ip_total_id <- paste('Ip', as.character(i), '_tatal_counts', sep="")
  input_id <- paste('Input', as.character(i), sep="")
  input_total_id <- paste('Input', as.character(i), '_tatal_counts', sep="")
  trt_m6a_ratio <- paste('m6a_ratio_trt_', as.character(i), sep="")
  
  trt_ratio <- (as.numeric(treat_data[ , ip_id])/as.numeric(treat_data[ , ip_total_id]))/ (as.numeric(treat_data[ , input_id])/as.numeric(treat_data[ , input_total_id]))
  
  temp_col<- append(colnames(combined_trt_ratio_data), trt_m6a_ratio)
  combined_trt_ratio_data <- cbind(combined_trt_ratio_data, trt_ratio)
  colnames(combined_trt_ratio_data) <- temp_col
}

# write.table(combined_trt_ratio_data, file="Treat_combined_counts_ratio_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#########################################################################################################
###    Control_Treatment_Integrated  processing                                          
#########################################################################################################

combined_ctr_trt_ratio_data <- merge(combined_ctr_ratio_data, combined_trt_ratio_data, by.x="Id", by.y="Id", incomparables = NA)

combined_sequence_data <- combined_ctr_trt_ratio_data[ , c("Id", "sequence.x")]

write.table(combined_sequence_data, file="combined_sequence_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

combined_colnames <- colnames(combined_ctr_trt_ratio_data)

ratio_cols <- grep('Id|m6a_ratio_ctr_[0-9]*$|m6a_ratio_trt_[0-9]*$', combined_colnames)

combined_ctr_trt_ratio_data_integrated <- combined_ctr_trt_ratio_data[ , ratio_cols]

write.table(combined_ctr_trt_ratio_data_integrated, file="Control_Treat_combined_counts_ratio_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)







