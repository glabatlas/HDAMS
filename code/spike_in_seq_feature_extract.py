# *coding:UTF-8
import math
import os
###############################
##  One-Hot encoding
###############################
def one_hot(raw_seq):
    my_seq=str(raw_seq)
    temp_ls=[]
    for i in range(0,len(my_seq)):
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
##  Chemical encoding
#########################################
def chemical(raw_seq):
    my_seq=str(raw_seq)
    temp_ls=[]
    for i in range(0,len(my_seq)):
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
##  Cumulative frequency of nucleotides
#########################################
def cumulative_nucle(raw_seq, step):
    raw_seq=str(raw_seq).upper()
    A_count=[]
    T_count=[]
    C_count=[]
    G_count=[]
    for i in range (len(raw_seq)):
        if i%step==0:
            A_count.append(raw_seq.count('A', 0, i+1)/float(i+1))
            T_count.append(raw_seq.count('T', 0, i+1)/float(i+1))
            C_count.append(raw_seq.count('C', 0, i+1)/float(i+1))
            G_count.append(raw_seq.count('G', 0, i+1)/float(i+1))
    code=A_count+T_count+C_count+G_count
    return code

#########################################
##  Impurity encoding
#########################################
def purity(raw_seq, wd_size):
    raw_seq=str(raw_seq).upper()
    A_all_seq=raw_seq.count('A')/float(len(raw_seq))
    T_all_seq=raw_seq.count('T')/float(len(raw_seq))
    C_all_seq=raw_seq.count('C')/float(len(raw_seq))
    G_all_seq=raw_seq.count('G')/float(len(raw_seq))
    A_sub_seq=[]
    T_sub_seq=[]
    C_sub_seq=[]
    G_sub_seq=[]
    for i in range(0,len(raw_seq)-wd_size+1):
        sub_seq=raw_seq[i:i+wd_size]
        A_sub_seq.append(sub_seq.count('A')/float(wd_size))
        T_sub_seq.append(sub_seq.count('T')/float(wd_size))
        C_sub_seq.append(sub_seq.count('C')/float(wd_size))
        G_sub_seq.append(sub_seq.count('G')/float(wd_size))

    clarity_feature=[]
    for k in range(len(A_sub_seq)):
        kl_score=A_sub_seq[k]*math.log(A_sub_seq[k]/(A_all_seq+1)+0.01, 2)+\
                 T_sub_seq[k]*math.log(T_sub_seq[k]/(T_all_seq+1)+0.01, 2)+\
                 C_sub_seq[k]*math.log(C_sub_seq[k]/(C_all_seq+1)+0.01, 2)+\
                 G_sub_seq[k]*math.log(G_sub_seq[k]/(G_all_seq+1)+0.01, 2)

        clarity_feature.append(kl_score)
    return clarity_feature

#########################################
##  Sequence-entropy-info encoding
#########################################
def entropy_inf(sub_seq):
    sub_seq=str(sub_seq).upper()
    A_sub_seq=sub_seq.count('A')/float(len(sub_seq))
    T_sub_seq=sub_seq.count('T')/float(len(sub_seq))
    C_sub_seq=sub_seq.count('C')/float(len(sub_seq))
    G_sub_seq=sub_seq.count('G')/float(len(sub_seq))

    temp_list=[A_sub_seq,T_sub_seq,C_sub_seq,G_sub_seq]
    entro_inf=0
    for x in temp_list:
        entro_inf+= -(x*math.log(x+0.01,2))
    return entro_inf

#########################################
##  Sequence-entropy-energy encoding
#########################################
def entropy_energy(raw_seq, wd_size):
    if len(raw_seq)<3*wd_size:
        return []
    else:
        center_pos=len(raw_seq)//2
        flank=(wd_size)//2
        center_seq=raw_seq[center_pos-flank:center_pos+flank+1]
        left_sub_seq=raw_seq[0:center_pos-flank]
        right_sub_temp=list(raw_seq[center_pos+flank+1:])
        right_sub_temp.reverse()

        right_sub_seq="".join(right_sub_temp)
        center_seq_entro= entropy_inf(center_seq)

        left_feature=[]
        right_feature=[]
        for i in range(len(left_sub_seq)-wd_size+1):
            slide_seq=left_sub_seq[i:i+wd_size]
            dist=len(left_sub_seq)-wd_size-i+1

            ifeature_lf=entropy_inf(slide_seq)
            join_seq=slide_seq+center_seq
            join_feature=entropy_inf(join_seq)

            entropy_dis_lf=(center_seq_entro + ifeature_lf - join_feature)/math.sqrt(dist)
            left_feature.append(entropy_dis_lf)

        for j in range(len(right_sub_seq)-wd_size+1):
            slide_seq=right_sub_seq[j:j+wd_size]

            dist=len(right_sub_seq)-wd_size-j+1
            ifeature_rg=entropy_inf(slide_seq)
            join_seq=slide_seq+center_seq
            join_feature=entropy_inf(join_seq)

            entropy_dis_rg=(center_seq_entro + ifeature_rg - join_feature)/math.sqrt(dist)
            right_feature.append(entropy_dis_rg)

        comb_feature=left_feature+right_feature
        return comb_feature

################################
########      main       #######
################################
if __name__ == "__main__":
    LEN=151   #define the length of coding sequence
    fread_seq=open('./training_data/design_oneline.txt', 'r')
    if not os.path.exists('./temp'):
        os.makedirs('./temp')

    fout_seq_1=open('./temp/seq_feature_one_hot.txt', 'w')
    fout_seq_2=open('./temp/seq_feature_chemical.txt', 'w')
    fout_seq_3=open('./temp/seq_feature_cumulative_nuclein.txt', 'w')
    fout_seq_4=open('./temp/seq_feature_clarity.txt', 'w')
    fout_seq_5=open('./temp/seq_feature_entropy_energy.txt', 'w')
    fout_seq_6=open('./temp/seq_feature_combined.txt', 'w')

    fread_seq_label=open('./training_data/seq_beta_avg_label_new.txt', 'r')
    fread_seq_label.readline()
    dic_lable={}
    for line in fread_seq_label:
        temp=str(line).strip().split('\t')
        key=temp[0]
        value=temp[1:]
        dic_lable.setdefault(key, value)
    fread_seq_label.close()

    dic_comb={}
    dic_seq={}
    for line in fread_seq:
        seq=line.strip()
        if seq[0]=='>':
            key=seq[1:]
        else:
            str_seq=seq[4:]
            seq_length=len(str_seq)
            if seq_length>=LEN:
                ancor_a=int((seq_length-LEN)/2)-1
                ancor_b=ancor_a+LEN
                str_seq=str_seq[ancor_a:ancor_b]
                dic_seq.setdefault(key, str_seq)
    fread_seq.close()

#######################################
########    One-hot 604           ####
#######################################
    dic_one_hot={}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        one_hot_code=one_hot(m6a_seq)
        dic_one_hot.setdefault(seq_id,one_hot_code)
        dic_comb.setdefault(seq_id,one_hot_code)
    print ('the one-hot coding features is:\t%d'%len(one_hot_code))

    for key in sorted(dic_lable.keys()):
        ls_value= dic_one_hot[key]
        seq_label=dic_lable[key][1]
        fout_seq_1.writelines('%s\t%s'%(seq_label, key))
        for vx in ls_value:
            fout_seq_1.writelines('\t%d'%vx)
        fout_seq_1.writelines('\n')
    fout_seq_1.close()

#######################################
########    chemical code   453    ####
#######################################
    dic_che_code={}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        chemical_code=chemical(m6a_seq)
        dic_che_code.setdefault(seq_id,chemical_code)
        ifeature=dic_comb[seq_id]
        ifeature.extend(chemical_code)
        dic_comb.setdefault(seq_id, ifeature)
    print ('the chemical coding features is:\t%d'%len(chemical_code))

    for key in sorted(dic_lable.keys()):
        ls_value= dic_che_code[key]
        seq_label=dic_lable[key][1]
        fout_seq_2.writelines('%s\t%s'%(seq_label,key))
        for vx in ls_value:
            fout_seq_2.writelines('\t%d'%vx)
        fout_seq_2.writelines('\n')
    fout_seq_2.close()

########################################
########    nuclein code    204     ####
########################################
    dic_nucle_code={}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        nucle_code=cumulative_nucle(m6a_seq, 3)
        dic_nucle_code.setdefault(seq_id,nucle_code)
        ifeature=dic_comb[seq_id]
        ifeature.extend(nucle_code)
        dic_comb.setdefault(seq_id, ifeature)
    print ('the nuclein coding features is:\t%d'%len(nucle_code))
    for key in sorted(dic_lable.keys()):
        ls_value= dic_nucle_code[key]
        seq_label=dic_lable[key][1]
        fout_seq_3.writelines('%s\t%s'%(seq_label,key))
        for vx in ls_value:
            fout_seq_3.writelines('\t%.4f'%vx)
        fout_seq_3.writelines('\n')
    fout_seq_3.close()

#######################################
########     impurity code  147   #####
#######################################
    dic_clarity_code={}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        clarity_code=purity(m6a_seq, 5)
        dic_clarity_code.setdefault(seq_id,clarity_code)
        # print len(clarity_code)
        ifeature=dic_comb[seq_id]
        ifeature.extend(clarity_code)
        dic_comb.setdefault(seq_id, ifeature)
    print ('the clarity coding features is:\t%d'%len(clarity_code))

    for key in sorted(dic_lable.keys()):
        ls_value= dic_clarity_code[key]
        seq_label=dic_lable[key][1]
        fout_seq_4.writelines('%s\t%s'%(seq_label, key))

        for vx in ls_value:
            fout_seq_4.writelines('\t%.4f'%vx)
        fout_seq_4.writelines('\n')
    fout_seq_4.close()

#######################################
####    entropy_energy code  108   ####
#######################################
    dic_energy_code={}
    for seq_id in dic_seq:
        m6a_seq=dic_seq[seq_id]
        energy_code=entropy_energy(m6a_seq, 15)
        dic_energy_code.setdefault(seq_id, energy_code)

        ifeature=dic_comb[seq_id]
        ifeature.extend(energy_code)
        dic_comb.setdefault(seq_id, ifeature)
        # print len(dic_comb[seq_id])
    print ('the energy coding features is:\t%d'%len(energy_code))

    for key in sorted(dic_lable.keys()):
        ls_value= dic_energy_code[key]
        seq_label=dic_lable[key][1]
        fout_seq_5.writelines('%s\t%s'%(seq_label,key))

        for vx in ls_value:
            fout_seq_5.writelines('\t%.4f'%vx)
        fout_seq_5.writelines('\n')
    fout_seq_5.close()

#################################
#####   combine_code         ####
#################################
    for key in sorted(dic_lable.keys()):
        ls_value= dic_comb[key]
        seq_label=dic_lable[key][1]
        fout_seq_6.writelines('%s\t%s'%(seq_label, key))
        for vx in ls_value:
            fout_seq_6.writelines('\t%.3f'%vx)
        fout_seq_6.writelines('\n')
    print ('the final combined coding features is:\t%d'%len(ls_value))
    fout_seq_6.close()
