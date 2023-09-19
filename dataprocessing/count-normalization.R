#!/sur/bin/Rscript
setwd("./HDAMS")

args <- commandArgs()
ctr_file <-args[6]
contr_group_num <-as.integer(args[7])

trt_file <-args[8]
treat_group_num <-as.integer(args[9])

exomepeak_sites_volap25  <- args[10]



contr_data <-read.table(ctr_file, header = TRUE)

treat_data <-read.table(trt_file, header = TRUE)

peak_data <-read.table(exomepeak_sites_volap25, header=FALSE)

event_id_ctr <-paste(contr_data[ , "peakname"], contr_data[ , "chrom"], contr_data[ , "strand"], contr_data[ , "start"], contr_data[ , "end"], sep="_")

event_id_trt <-paste(treat_data[ , "peakname"], treat_data[ , "chrom"], treat_data[ , "strand"], treat_data[ , "start"], treat_data[ , "end"], sep="_")

event_id_peak <-peak_data[ , 1]

contr_data <-cbind(Id=event_id_ctr, contr_data)
contr_peak_index <-match(event_id_peak, contr_data[,"Id"], nomatch = FALSE)
contr_data <-contr_data[contr_peak_index, ]


treat_data <-cbind(Id=event_id_trt, treat_data)
treat_peak_index <-match(event_id_peak, treat_data[,"Id"], nomatch = FALSE)
treat_data <-treat_data[treat_peak_index,]

combined_data <-merge(contr_data, treat_data, by.x="Id", by.y="Id", incomparables = NA)

###############################################################################################################
###  Input count processing 
###############################################################################################################

combined_data_cols <-colnames(combined_data)

input_cols <- combined_data_cols[grep('Id|Input[0-9]*[^_]*$',combined_data_cols)]

combined_input_counts_raw <- combined_data[ ,input_cols]

S1_row <- dim(combined_input_counts_raw)[1]
K1_col <- dim(combined_input_counts_raw)[2]

flag_1 <- rep(TRUE, times=S1_row)

for (i in seq(1:S1_row))
{   
  if (sum(c(combined_input_counts_raw[i, 2:K1_col])<1) >0 | mean(as.numeric(combined_input_counts_raw[i,2:K1_col]))<20)
  { 
    flag_1[i]<-FALSE
  }
  
}

combined_input_counts_raw <- combined_input_counts_raw[flag_1, ]

combined_input_counts <- combined_input_counts_raw

S1_row <- dim(combined_input_counts)[1]
K1_col <- dim(combined_input_counts)[2]

# write.table(combined_input_counts, file="siControl_Treatment_combined_input_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Norm_factor <-rep(0, times=S1_row)

for (i in seq(1:S1_row))
{ 
  r_fac=1
  
  for(j in seq(from=2, to=K1_col, by=1))
  {
    if (combined_input_counts[i,j]!=0)
      
    { r_fac <- r_fac*combined_input_counts[i,j] }
  }
  
  eff_len <- sum(c(combined_input_counts[i, 2:K1_col])!=0)
  
  Norm_factor[i] <- r_fac^(1/eff_len)
}

for (i in seq(1:S1_row))
{
  combined_input_counts[i, 2:K1_col] <-as.numeric(combined_input_counts[i, 2:K1_col])/Norm_factor[i] 
}

for (j in seq(from=2, to=K1_col, by=1))
{
  median_factor <- median(combined_input_counts[ ,j])
  # cat(median_factor)
  # cat("\n")
  combined_input_counts_raw[ ,j] <- floor(combined_input_counts_raw[ ,j]*median_factor)
}

write.table(combined_input_counts_raw, file="siControl_Treatment_combined_input_counts_normalized.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#############################################
####  IP counts normalization   #############
#############################################

combined_data_cols <-colnames(combined_data)

ip_cols <- combined_data_cols[grep('Id|Ip[0-9]*[^_]*$' ,combined_data_cols)]

combined_ip_counts_raw <- combined_data[ ,ip_cols]

S2_row <-dim(combined_ip_counts_raw)[1]
K2_col <-dim(combined_ip_counts_raw)[2]

flag_2<-rep(TRUE, times=S2_row)

for (i in seq(1:S2_row))
{   
  if (sum(c(combined_ip_counts_raw[i,2:K2_col])<1) >0 | mean(as.numeric(combined_ip_counts_raw[i, 2:K2_col])) <20)
  { 
    flag_2[i]<-FALSE
  }
}

combined_ip_counts_raw <- combined_ip_counts_raw[flag_2, ]

combined_ip_counts <- combined_ip_counts_raw

S2_row <-dim(combined_ip_counts)[1]
K2_col <-dim(combined_ip_counts)[2]

# write.table(combined_ip_counts, file="siControl_Treatment_combined_ip_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Norm_factor <-rep(0, times=S2_row)

for (i in seq(1:S2_row))
{
  r_fac=1
  
  for(j in seq(from=2, to=K2_col, by=1))
  {
    if (combined_ip_counts[i,j]!=0)
    
      { r_fac <- r_fac*combined_ip_counts[i,j] }
  }
  
  eff_len <-sum(c(combined_ip_counts[i, 2:K2_col])!=0)
  
  Norm_factor[i] <-r_fac^(1/eff_len)
}

for (i in seq(1:S2_row))
{
  combined_ip_counts[i, 2:K2_col] <- as.numeric(combined_ip_counts[i, 2:K2_col])/Norm_factor[i]
  
}

for (j in seq(from=2, to=K2_col, by=1))
{
  median_factor <- median(combined_ip_counts[ ,j])
  # cat(median_factor)
  # cat("\n")
  
  combined_ip_counts_raw[ ,j] <-floor(combined_ip_counts_raw[ ,j]*median_factor)
  
}

write.table(combined_ip_counts_raw, file="siControl_Treatment_combined_ip_counts_normalized.txt", sep = "\t", quote = FALSE, row.names = FALSE)


######################################################################
##     Integrated counts formatted 
######################################################################

ip_cols_normal <- colnames(combined_ip_counts_raw)
input_cols_normal <- colnames(combined_input_counts_raw)

for ( m in seq(from=2, to=length(ip_cols_normal), by=1))
{
  if (m <=contr_group_num+1)
   
    { ip_cols_normal[m] <- paste("Ip",as.character(m-1), "_ctr", sep="")}
  
  else
   
    { ip_cols_normal[m] <- paste("Ip",as.character(m-contr_group_num-1), "_trt", sep="") }
}

for ( n in seq(from=2, to=length(input_cols_normal), by=1))
{
  if (n <=contr_group_num+1)
  
    { input_cols_normal[n] <- paste("Input", as.character(n-1),"_ctr", sep="")}
 
  else
  
    { input_cols_normal[n] <- paste("Input", as.character(n-contr_group_num-1), "_trt", sep="")}
}

colnames(combined_ip_counts_raw) <-ip_cols_normal
colnames(combined_input_counts_raw) <-input_cols_normal

combined_normalized_count_data <- merge(combined_ip_counts_raw, combined_input_counts_raw, by.x="Id", by.y="Id", incomparables = NA)


write.table(combined_normalized_count_data, file="siControl_Teatment_combined_normalized_count_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)




