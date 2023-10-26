import re
import os
#####################
def transcript_file_read(filename):
    '''
    input:
    Ensemble_transcript_name   chrom   strand   GeneStart   GeneEnd CdsStart    CdsEnd  exonCount   exonStart   exonEnd socre   Ensemble_genename  type    genename
    ENST00000473358 chr1    +       29553   31097   31097   31097   3       29553,30563,30975,      30039,30667,31097, 0       ENSG00000243485 none    none    -1,-1,-1,
    ENST00000469289 chr1    +       30266   31109   31109   31109   2       30266,30975,    30667,31109,    0 ENSG00000243485  none    none    -1,-1,

    output: dict[Ensemble_gene_id]=[exonstart,exonend,Ensemble_trans_id,chrom,strand]
    '''
    dictout = {}
    f = open(filename,'r')
    for str_x in f:
        str_x = str_x.strip("\n")
        list_x = str_x.split('\t')
        if str_x[0]=='#' or list_x[0]=="Ensemble_transcript_name":
            #print(list_x)
            index_1 = list_x.index("Ensemble_transcript_name")
            index_2 = list_x.index("chrom")
            index_3 = list_x.index("strand")
            index_4 = list_x.index("GeneStart")
            index_5 = list_x.index("GeneEnd")
            index_6 = list_x.index("CdsStart")
            index_7 = list_x.index("CdsEnd")
            index_8 = list_x.index("exonCount")
            index_9 = list_x.index("exonStart")
            index_10 = list_x.index("exonEnd")
            index_11 = list_x.index("Ensemble_gene_name")
            continue
        Ensemble_trans_id = list_x[index_1]
        chrom = list_x[index_2]
        strand = list_x[index_3]
        genestart = list_x[index_4]
        geneend = list_x[index_5]
        cdsstart = list_x[index_6]
        cdsend = list_x[index_7]
        startlist = list_x[index_9].split(',')[:-1]
        endlist = list_x[index_10].split(',')[:-1]
        Ensemble_gene_id = list_x[index_11]
        for i in range(len(startlist)):
            try:
                dictout[Ensemble_gene_id].append([int(startlist[i]),int(endlist[i]),Ensemble_trans_id,chrom,strand,cdsstart,cdsend,genestart,geneend])
                #dictout[Ensemble_gene_id]=[[int(startlist[i]),int(endlist[i]),Ensemble_trans_id,chrom,strand]]
            except:
                dictout[Ensemble_gene_id]=[[int(startlist[i]),int(endlist[i]),Ensemble_trans_id,chrom,strand,cdsstart,cdsend,genestart,geneend]]
                #dictout[Ensemble_gene_id].append([int(startlist[i]),int(endlist[i]),Ensemble_trans_id,chrom,strand])
    return dictout

####################################
def same_gene_transcript_merge(dictinput,result_filename,chrom2read_dict,path):
    '''
    output: dict[Ensemble_gene_id]=[exonstart,exonend,Ensemble_trans_id,chrom,strand]
    '''
    writefile = open(path+result_filename,'a')
    for keys in dictinput.keys():
        chrom = dictinput[keys][0][3]
        strand = dictinput[keys][0][4]
        dictinput[keys].sort(key = lambda x:x[0])
        exonlist = dictinput[keys]
        min_cds = min([int(x[5]) for x in exonlist])
        max_cds = max([int(x[6]) for x in exonlist])
        min_gene = min([int(x[7]) for x in exonlist])
        max_gene = max([int(x[8]) for x in exonlist])
        strand_plus = '\t'.join([str(min_gene),str(min_cds),str(max_cds),str(max_gene)])
        strand_minus = '\t'.join([str(max_gene),str(max_cds),str(min_cds),str(min_gene)])
        #print(exonlist)
        i = 0
        while (i < len(exonlist)-1):
            if exonlist[i][1] > exonlist[i+1][0] and exonlist[i+1][1] > exonlist[i][0]:
                exonlist[i][1] = max(exonlist[i][1],exonlist[i+1][1])
                exonlist[i][0] = min(exonlist[i][0],exonlist[i+1][0])
                exonlist.pop(i+1)
            else:
                i+=1
        endlist = []
        startlist = []
        for j in range(len(exonlist)):
            endlist.append(str(exonlist[j][1]))
            startlist.append(str(exonlist[j][0]))
        if strand=="+":
            StrToWrite = keys+'\t'+chrom+'\t'+','.join(startlist)+','+'\t'+','.join(endlist)+','+'\t'+strand+'\t'+strand_plus+'\n'
        else:
            StrToWrite = keys+'\t'+chrom+'\t'+','.join(startlist)+','+'\t'+','.join(endlist)+','+'\t'+strand+'\t'+strand_minus+'\n'
        DRACH_sub_result_list = sequecne_stract(StrToWrite=StrToWrite,chrom2read_dict=chrom2read_dict)
        bed_format_resultwrite(writefile=writefile,DRACH_sub_result_list=DRACH_sub_result_list)


def bed_format_resultwrite(writefile,DRACH_sub_result_list):
    for str_x in DRACH_sub_result_list:
        writefile.write(str_x)

def sequecne_stract(StrToWrite,chrom2read_dict):
    """
    ENSG00000157191 chr1    16767166,16767660,16770126,16774364,16774554,16774743,16775587,16782015,16785336, 16767348,16768114,16770227,16774469,16774636,16774845,16778510,16782388,16786573,        +       16767166  16767256 16786572        16786573
    ENSG00000169504 chr1    25071847,25119667,25124232,25140584,25153500,25166350,25167263, 25072116,25119704,25124454,25140710,25153607,25166532,25170815,    +       25071847        25072044        25167428        25170815
    """
    str_x = StrToWrite.strip("\n")
    list_x = str_x.split("\t")
    ###########
    out_result_list = []
    ###########
    Ensemble_gene_id = list_x[0]
    chrom = list_x[1]
    exon_start_list = [int(x) for x in list_x[2].split(",")[:-1]]
    exon_end_list = [int(x) for x in list_x[3].split(",")[:-1]]
    strand = list_x[4]
    gene_start = list_x[5]
    cds_start = list_x[6]
    cds_end = list_x[7]
    gene_end = list_x[8]
    ###########
    for exon_num_i in range(len(exon_start_list)):
        exon_start_i = int(exon_start_list[exon_num_i])
        exon_end_i = int(exon_end_list[exon_num_i])
        out_result_list_i = DRACT_sub_sequecne_abstract(
            chrom2read_dict=chrom2read_dict,
            chrom=chrom,
            start=exon_start_i,
            end=exon_end_i,
            strand=strand,
            ensemble_gene_id=Ensemble_gene_id,
            exon_id=exon_num_i
            )
            #####################################
        out_result_list += out_result_list_i
        ###########
    return out_result_list

def chrom_to_read_dict_make(filename):
    f = open(filename,'r')
    outdict = {}
    for str_x in f:
        str_x = str_x.strip("\n")
        #list_x = str_x.split("\t")
        #outdict[list_x[0]] = list_x[1]
        if str_x[0] == ">":
            chrom = str_x[1:]
            outdict[chrom] = []
        else:
            outdict[chrom].append(str_x)
    oneline_dict = {}
    for chrom in outdict.keys():
        onelineread = "".join(outdict[chrom])
        oneline_dict[chrom] = onelineread
    return oneline_dict

def DRACT_sub_sequecne_abstract(chrom2read_dict,chrom,start,end,strand,ensemble_gene_id,exon_id):
    out_result_list = []
    try:
        tmpread = chrom2read_dict[chrom][start:end]
    except:
        return out_result_list
    #####################
    for i in range(len(tmpread) - 151):
        tmp_read_i = tmpread[i:i+151]
        start_i = start+i
        end_i = start_i+151
        ######################
        sitename = "_".join([ensemble_gene_id,"exon",str(exon_id+1),str(i+1)])
        ######################
        if strand == "+":
            tmp_read_convert = tmp_read_i
        else:
            tmp_read_convert = reverse_make(read=tmp_read_i)
        match_DRACH = motif_parten_match(sequence=tmp_read_convert[75-2:75+3],motif_parten="[GAT][GA]AC[ATC]")
        if match_DRACH == "True":
            list2write = [chrom,str(start_i),str(end_i),sitename,ensemble_gene_id,strand,tmp_read_convert]
            str2write = "\t".join(list2write) + "\n"
            out_result_list.append(str2write)
        else:
            pass
    return out_result_list

def motif_parten_match(sequence,motif_parten):
    motif_list = re.findall(motif_parten,sequence)
    if motif_list == []:
        return "False"
    else:
        return "True"

def reverse_make(read):
    read.upper()
    convert_dict = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
    ####################
    tmplist = []
    for i in range(len(read)):
        base = read[i]
        convert_base = convert_dict[base]
        tmplist.append(convert_base)
    tmplist.reverse()
    tmpstr = "".join(tmplist)
    return tmpstr

###########################################################

def multicov_pipeline(input_bam_list,ip_bam_list,bed,path):
    input_bam_str = " ".join(input_bam_list)
    ip_bam_str = " ".join(ip_bam_list)
    str2ext = " ".join(["bedtools multicov -q 20 -bams",input_bam_str,ip_bam_str,"-bed",bed,">>",path+"tmp_multicov_result.txt"])
    os.system(str2ext)

def samtools_pipeline(input_bam_list,ip_bam_list,path):
    input_count_list = []
    for bam in input_bam_list:
        str2ext = " ".join(["samtools view -q 20",bam,"| wc -l >>",path+"tmp_input_bam_total_read.txt"])
        print(str2ext)
        os.system(str2ext)
    for bam in ip_bam_list:
        str2ext = " ".join(["samtools view -q 20",bam,"| wc -l >>",path+"tmp_ip_bam_total_read.txt"])
        print(str2ext)
        os.system(str2ext)
    ########################
    input_count_list = mapped_depth_count(filename=path+'tmp_input_bam_total_read.txt')
    ip_count_list = mapped_depth_count(filename=path+'tmp_ip_bam_total_read.txt')
    ########################
    return input_count_list,ip_count_list

def mapped_depth_count(filename):
    f = open(filename,'r')
    countlist = []
    for str_x in f:
        str_x = str_x.strip("\n")
        count = str_x
        countlist.append(count)
    return countlist

def resultwrite(multicov_result,input_count_list,ip_count_list,out_count_result,path):
    f = open(multicov_result,'r')
    d = open(path+out_count_result,'a')
    for str_x in f:
        str_x = str_x.strip("\n")
        outstr = "\t".join([str_x]+input_count_list+ip_count_list) + "\n"
        d.write(outstr)

def args_make():
    import argparse
    parser = argparse.ArgumentParser(description='DRACH motif sequence abstract')
    parser.add_argument('--transcript_filename', required=True, help="transcript file name")
    parser.add_argument('--fasta_filename', required=True, help="fasta file name")
    parser.add_argument('--bedformat_result_filename', required=True, help="bedformat result file name")
    parser.add_argument('--input_bams', required=True, help="input bams:bam1,bam2,bam3")
    parser.add_argument('--ip_bams', required=True, help="ip bams:bam1,bam2,bam3")
    parser.add_argument('--count_result_filename', required=True, help="result file name")
    parser.add_argument('--tmp_path', required=True, help="result file name")
    args = parser.parse_args()
    return args

if __name__=='__main__':
    import sys
    args = args_make()
    transcript_filename = args.transcript_filename
    fasta_filename = args.fasta_filename
    bedformat_result_filename = args.bedformat_result_filename
    input_bams = args.input_bams
    ip_bams = args.ip_bams
    count_result_filename = args.count_result_filename
    tmp_path = args.tmp_path
    ####################################
    ####################################
    transcript_list = transcript_file_read(filename=transcript_filename)
    chrom2read_dict = chrom_to_read_dict_make(filename=fasta_filename)
    same_gene_transcript_merge(
        dictinput=transcript_list,
        result_filename=bedformat_result_filename,
        chrom2read_dict=chrom2read_dict,
        path=tmp_path
        )
    print("DRACH motif sequence generated form transcript is abstracted!!!")
    ####################################
    input_bam_list = input_bams.split(",")
    ip_bam_list = ip_bams.split(",")
    multicov_pipeline(
        input_bam_list=input_bam_list,
        ip_bam_list=ip_bam_list,
        bed=bedformat_result_filename,
        path=tmp_path
        )
    print("multi coverage is calculated!!!")
    ####################################
    input_count_list,ip_count_list = samtools_pipeline(
        input_bam_list=input_bam_list,
        ip_bam_list=ip_bam_list,
        path=tmp_path
        )
    ####################################
    resultwrite(
        multicov_result="tmp_multicov_result.txt",
        input_count_list=input_count_list,
        ip_count_list=ip_count_list,
        out_count_result=count_result_filename,
        path=tmp_path
        )
    print("total coverage was added!!!")
    