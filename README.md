# Summary 

HDAMS is a hierarchical statistic model to identify differential m6A sites among sample groups by considering the effect of antibody-specificty of sequence. 

# Dependencies 

  python 2.7;  numpy 1.16.6;  scipy 1.2.1;  pandas 0.24.2;  scikit-learn 0.20.3.

# Installation 
git this project directory

cd ./HDAMS/

# Usage
## Data preprocessing 
* **Obtain original m6a ratio:** original count data of all m6a site-candidate sequences (with DRACH motif) extracted from BAM file will be used to calculate the original m6a ratio of each sequence.
  
        >> Rscript ./dataprocessing/Origin-ratio-calculation.R control_group_file.txt num_of_smaples_in_control treat_group_file.txt num_of_smaples_in_treat
  
* **Count data normalization:** to obtain the normalized IP and INPUT counts. The sites: (1) seqeucnces overlap with all exomepeak calling peaks < 25bp; (2)average IP or INPUT counts < 20; will be automatically filtered out before normalization. 
  
        >> Rscript ./dataprocessing/count-normalization.R control_group_file.txt num_of_smaples_in_control treat_group_file.txt num_of_smaples_in_treat selected_seqs_list_overlap_with_exomepeak_calling.txt

* **Obtain differential analysis format:** differential format data will be used to input HDAMS differential model to perform differential analysis of sites.

        >> Rscript ./dataprocessing/get-diff-format-data.R control_treat_combined_normalized_count_file.txt control_treat_combined_m6a_ratio_data.txt predicted_antibody_specificity_level.txt
  
## Site-specific antibody specificity prediction
* **Predict site-specific antibody specificity:** input the sequences (center-position is the condidate site, with 151bp) and multiple types features are automatically extracted to predict its specificty.  

      >> python ./code/predict_specificity.py --input input_sequence_file.txt --length seq_len --output output_file.txt  

* **Re-train the default ensemble prediction model:** if you want to re-train the default fitted model, we also provide the training code and can be used directly as follows. (Note: new fitted model will be used and may have minor inconsistent comparing with the default model, as randomly tree-based feature selection strategy was used.)

        >> python ./code/model-fitting.py   

## Differential analysis
Differential analysis based on the differential data format handled by data processing. It can be performed by,

    >> python ./code/hdams-diff.py --input input_diff_format_data.txt --fdr threshold --output output_results.txt


**Note:** (1) to predict your sequences' antibody specificities, the sequence data file should put into the fold of ./predict_data/; (2) sequence should be with lenght of 151bp and center-position is the candidate site; (3) re-train prediction model will over-write the defalut model and will use the new model to predict the seqences' specificities. 

