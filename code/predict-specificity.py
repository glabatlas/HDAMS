import math
import os, sys, argparse
from sklearn.externals import joblib
# import joblib
import pandas as pd
import numpy as np
###############################
##  One-hot Encoding
###############################
def one_hot(raw_seq):
    my_seq=str(raw_seq)
    temp_ls=[]
    for i in range(0, len(my_seq)):
        if my_seq[i]=='A':
            temp_ls.extend([1, 0, 0, 0])
        elif my_seq[i]=='C':
            temp_ls.extend([0, 1, 0, 0])
        elif my_seq[i]=='G':
            temp_ls.extend([0, 0, 1, 0])
        elif my_seq[i]=='T':
            temp_ls.extend([0, 0, 0, 1])
    return temp_ls

#########################################
##  Chemical Encoding
#########################################
def chemical(raw_seq):
    my_seq=str(raw_seq)
    temp_ls=[]
    for i in range(0, len(my_seq)):
        if my_seq[i]=='A':
            temp_ls.extend([1, 1, 1])
        elif my_seq[i]=='C':
            temp_ls.extend([0, 1, 0])
        elif my_seq[i]=='G':
            temp_ls.extend([1, 0, 0])
        elif my_seq[i]=='T':
            temp_ls.extend([0, 0, 1])
    return temp_ls

#########################################
##  Cumulative Nuclein Encoding
#########################################
def cumulative_nucle(raw_seq, step):
    raw_seq=str(raw_seq).upper()
    A_count=[]
    T_count=[]
    C_count=[]
    G_count=[]
    for i in range(len(raw_seq)):
        if i%step==0:
            A_count.append(raw_seq.count('A', 0, i+1)/float(i+1))
            T_count.append(raw_seq.count('T', 0, i+1)/float(i+1))
            C_count.append(raw_seq.count('C', 0, i+1)/float(i+1))
            G_count.append(raw_seq.count('G', 0, i+1)/float(i+1))
    code=A_count+T_count+C_count+G_count
    return  code

#########################################
##  Purity Encoding
#########################################
def purity(raw_seq, wd_size):
    raw_seq=str(raw_seq).upper()
    A_all_seq= raw_seq.count('A')/float(len(raw_seq))
    T_all_seq= raw_seq.count('T')/float(len(raw_seq))
    C_all_seq= raw_seq.count('C')/float(len(raw_seq))
    G_all_seq= raw_seq.count('G')/float(len(raw_seq))
    A_sub_seq=[]
    T_sub_seq=[]
    C_sub_seq=[]
    G_sub_seq=[]
    for i in range(0, len(raw_seq)-wd_size+1):
        sub_seq= raw_seq[i:i+wd_size]
        A_sub_seq.append(sub_seq.count('A')/float(wd_size))
        T_sub_seq.append(sub_seq.count('T')/float(wd_size))
        C_sub_seq.append(sub_seq.count('C')/float(wd_size))
        G_sub_seq.append(sub_seq.count('G')/float(wd_size))

    clarity_feature=[]
    for k in range(len(A_sub_seq)):
        kl_score= A_sub_seq[k]*math.log(A_sub_seq[k]/(A_all_seq+1)+0.01, 2)+\
                  T_sub_seq[k]*math.log(T_sub_seq[k]/(T_all_seq+1)+0.01, 2)+\
                  C_sub_seq[k]*math.log(C_sub_seq[k]/(C_all_seq+1)+0.01, 2)+\
                  G_sub_seq[k]*math.log(G_sub_seq[k]/(G_all_seq+1)+0.01, 2)
        clarity_feature.append(kl_score)
    return clarity_feature

#########################################
##  Sequence-Entropy-Info Encoding
#########################################
def entropy_inf(sub_seq):
    sub_seq= str(sub_seq).upper()
    A_sub_seq= sub_seq.count('A')/float(len(sub_seq))
    T_sub_seq= sub_seq.count('T')/float(len(sub_seq))
    C_sub_seq= sub_seq.count('C')/float(len(sub_seq))
    G_sub_seq= sub_seq.count('G')/float(len(sub_seq))

    temp_list= [A_sub_seq, T_sub_seq, C_sub_seq, G_sub_seq]
    entro_inf= 0
    for x in temp_list:
        entro_inf+= -(x*math.log(x+0.01, 2))
    return entro_inf

#########################################
##  Sequence-Entropy-Energy Encoding
#########################################
def entropy_energy(raw_seq, wd_size):
    if len(raw_seq) < 3*wd_size:
        return []
    else:
        center_pos= len(raw_seq)//2
        flank= (wd_size)//2
        center_seq= raw_seq[center_pos-flank:center_pos+flank+1]
        left_sub_seq= raw_seq[0:center_pos-flank]
        right_sub_temp= list(raw_seq[center_pos+flank+1:])
        right_sub_temp.reverse()

        right_sub_seq= "".join(right_sub_temp)
        center_seq_entro= entropy_inf(center_seq)

        left_feature=[]
        right_feature=[]
        for i in range(len(left_sub_seq)-wd_size+1):
            slide_seq= left_sub_seq[i:i+wd_size]
            dist= len(left_sub_seq)-wd_size-i+1

            ifeature_lf= entropy_inf(slide_seq)
            join_seq= slide_seq+center_seq
            join_feature= entropy_inf(join_seq)

            entropy_dis_lf= (center_seq_entro + ifeature_lf - join_feature)/math.sqrt(dist)
            left_feature.append(entropy_dis_lf)

        for j in range(len(right_sub_seq)-wd_size+1):
            slide_seq= right_sub_seq[j:j+wd_size]
            dist= len(right_sub_seq)-wd_size-j+1
            ifeature_rg= entropy_inf(slide_seq)
            join_seq= slide_seq+center_seq
            join_feature= entropy_inf(join_seq)

            entropy_dis_rg=(center_seq_entro + ifeature_rg - join_feature)/math.sqrt(dist)
            right_feature.append(entropy_dis_rg)

        comb_feature=left_feature+right_feature
        return comb_feature


def pred_model_load(model_name):
    model_list = os.listdir('./model/')
    select_model = model_name
    my_model = []
    if len(model_list) > 0 :
        for i in range(len(model_list)):
            if model_list[i] in select_model:
                i_model = joblib.load('./model/'+model_list[i])
                my_model.append(i_model)
    return my_model


if __name__ == "__main__":
    LEN=151   ## defining the length of coding sequence

    parser=argparse.ArgumentParser(description='please provide your sequence file!')
    parser.add_argument('--input', required=True, help='specify the input data file')
    parser.add_argument('--length', required=False, default=151, help='specify sequence length')
    parser.add_argument('--output', required=True, help='specify the output file')

    args=parser.parse_args()
    input_file=args.input

    LEN=int(args.length)
    output_file=args.output

    if len([input_file, output_file])<2:
        print("please specify your input and output file!")
        sys.exit(1)

    fread_seq=open('./predict_data/'+input_file, 'r')

    fout_seq_1=open('./predict_data/seq_feature_one_hot_new.txt','w')
    fout_seq_2=open('./predict_data/seq_feature_chemical_new.txt','w')
    fout_seq_3=open('./predict_data/seq_feature_cumulative_nuclein_new.txt','w')
    fout_seq_4=open('./predict_data/seq_feature_clarity_new.txt','w')
    fout_seq_5=open('./predict_data/seq_feature_entropy_energy_new.txt','w')
    fout_seq_6=open('./predict_data/seq_feature_combined_new.txt','w')

    dic_comb={}
    dic_seq={}
    for line in fread_seq:
        seq=line.strip().split('\t')
        key=seq[0]
        str_seq=str.upper(seq[1])
        seq_length=len(str_seq)
        if seq_length<LEN:
            continue
        elif seq_length>= LEN:
            ancor_a= int((seq_length-LEN)/2)
            ancor_b= ancor_a+LEN
            str_seq= str_seq[ancor_a:ancor_b]
            center_position= LEN//2

            if str_seq[center_position]=='A':
               dic_seq.setdefault(key, str_seq)
    fread_seq.close()

    #######################################
    ##    One-hot 604
    #######################################
    dic_one_hot={}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        one_hot_code=one_hot(m6a_seq)
        dic_one_hot.setdefault(seq_id, one_hot_code)
        dic_comb.setdefault(seq_id, one_hot_code)
    print('one-hot coding features:\t%d'% len(one_hot_code))

    for key in sorted(dic_seq.keys()):
        ls_value=dic_one_hot[key]
        fout_seq_1.writelines('%s'% key)
        for vx in ls_value:
            fout_seq_1.writelines('\t%d'% vx)
        fout_seq_1.writelines('\n')
    fout_seq_1.close()

    #######################################
    ########    chemical code   453
    #######################################
    dic_che_code={}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        chemical_code=chemical(m6a_seq)
        dic_che_code.setdefault(seq_id, chemical_code)

        ifeature=dic_comb[seq_id]
        ifeature.extend(chemical_code)
        dic_comb.setdefault(seq_id, ifeature)
    print ('chemical coding features:\t%d'% len(chemical_code))

    for key in sorted(dic_seq.keys()):
        ls_value= dic_che_code[key]
        fout_seq_2.writelines('%s'% key)
        for vx in ls_value:
            fout_seq_2.writelines('\t%d'%vx)
        fout_seq_2.writelines('\n')
    fout_seq_2.close()

    ########################################
    ########    nuclein code    204     ####
    ########################################
    dic_nucle_code= {}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        nucle_code=cumulative_nucle(m6a_seq, 3)
        dic_nucle_code.setdefault(seq_id, nucle_code)

        ifeature=dic_comb[seq_id]
        ifeature.extend(nucle_code)
        dic_comb.setdefault(seq_id, ifeature)
    print ('nuclein coding features:\t%d'% len(nucle_code))

    for key in sorted(dic_seq.keys()):
        ls_value= dic_nucle_code[key]
        fout_seq_3.writelines('%s'% key)
        for vx in ls_value:
            fout_seq_3.writelines('\t%.4f'% vx)
        fout_seq_3.writelines('\n')
    fout_seq_3.close()

    #######################################
    ########     purity code  147   ######
    #######################################
    dic_purity_code= {}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        purity_code=purity(m6a_seq, 5)
        dic_purity_code.setdefault(seq_id, purity_code)

        ifeature=dic_comb[seq_id]
        ifeature.extend(purity_code)
        dic_comb.setdefault(seq_id, ifeature)
    print ('clarity coding features:\t%d'% len(purity_code))

    for key in sorted(dic_seq.keys()):
        ls_value=dic_purity_code[key]
        fout_seq_4.writelines('%s'% key)
        for vx in ls_value:
            fout_seq_4.writelines('\t%.4f'% vx)
        fout_seq_4.writelines('\n')
    fout_seq_4.close()

    #######################################
    ########    entropy_energy code
    #######################################
    dic_energy_code= {}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        energy_code=entropy_energy(m6a_seq, 15)
        dic_energy_code.setdefault(seq_id, energy_code)
        ifeature=dic_comb[seq_id]
        ifeature.extend(energy_code)
        dic_comb.setdefault(seq_id, ifeature)
    print ('energy coding features:\t%d'% len(energy_code))

    for key in sorted(dic_seq.keys()):
        ls_value= dic_energy_code[key]
        fout_seq_5.writelines('%s'% key)
        for vx in ls_value:
            fout_seq_5.writelines('\t%.4f'% vx)
        fout_seq_5.writelines('\n')
    fout_seq_5.close()

    #################################
    #####   combine_code
    #################################
    for key in sorted(dic_seq.keys()):
        ls_value= dic_comb[key]
        fout_seq_6.writelines('%s'% key)
        for vx in ls_value:
            fout_seq_6.writelines('\t%.3f'% vx)
        fout_seq_6.writelines('\n')
    print('final combined coding features:\t%d'% len(ls_value))
    fout_seq_6.close()

    ##############################################################
    print('Sequences encoding finished!!')
    ##############################################################
    ###  using trained model to predict sequence specificity
    ###############################################################

    processed_data = pd.read_csv('./predict_data/seq_feature_combined_new.txt', sep='\t', header=None)
    print("the dims of predicted data are:({0[0]}, {0[1]})".format(processed_data.shape))
    precessed_data_feature = np.asarray(processed_data.iloc[:, 1:])
    feature_id_input = open('./training_data/selected_feature_index.txt', 'r')
    indices_seq_id = np.asarray(feature_id_input.readline().strip().split('\t'), dtype=int)
    print ("num of selected features is:{0}".format(len(indices_seq_id)))
    selected_processed_data = precessed_data_feature[:, indices_seq_id]
    ###################################################################
    ##  model with: svm-lr
    ##################################################################
    model_name_1 = ['svm.joblib', 'lr.joblib']
    predict_model_1 = pred_model_load(model_name_1)

    res_predict_model_0 = predict_model_1[0].predict_proba(selected_processed_data)
    res_predict_model_1 = predict_model_1[1].predict_proba(selected_processed_data)

    predict_avg_prob_1 = (np.asmatrix(res_predict_model_0) + np.asmatrix(res_predict_model_1))/2
    SN_1 = processed_data.iloc[:, 0].shape[0]

    fout_predict_data_1 = open('./predict_data/'+output_file,'w')
    fout_predict_data_1.writelines('seq_id\tprob_beta_low\tprob_beta_high\n')
    for i in range(SN_1):
        fout_predict_data_1.writelines('%s'%processed_data.iloc[:, 0][i])
        ls_prob_1=np.asarray(predict_avg_prob_1[i, :])
        for x in ls_prob_1[0,:]:
            fout_predict_data_1.writelines('\t%.3f'% x)
        fout_predict_data_1.writelines('\n')
    fout_predict_data_1.close()
    print("Finished!")
