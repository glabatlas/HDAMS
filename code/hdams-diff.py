import sys
import warnings, scipy, operator, argparse
from sklearn import linear_model
from multiprocessing import Pool
import numpy as np
from scipy import stats
from scipy.optimize import fmin_l_bfgs_b
from math import log
# from scipy.special import gamma
warnings.filterwarnings('ignore')
#############################################################
##  likelihood function for the differential model
#############################################################
def likelihood_function(x, *args):
    beta_1=x[0]
    beta_2=x[1]
    theta_0=x[2]

    var_group=args[0]
    var_m6a_psi=args[1]
    var_ip_counts=args[2]
    all_counts=args[3]
    var_effi_x0=args[4]
    # var_input_counts=all_counts-var_ip_counts
    f_sum=0.0

    for i in range(len(var_m6a_psi)):
        f_sum+= -1*(var_ip_counts[i]*(beta_1*var_group[i]*var_effi_x0 + beta_2*var_m6a_psi[i] + theta_0)- \
                    all_counts[i]*log(1+np.exp(beta_1*var_group[i]*var_effi_x0 + beta_2*var_m6a_psi[i] + theta_0)))

    return(f_sum)

#######################################################
###  One-ordered function of derivative for likelihood
#######################################################
def likelihood_fun_der(x, *args):
    beta_1 = x[0]
    beta_2 = x[1]
    theta_0 = x[2]

    var_group = args[0]
    var_m6a_psi = args[1]
    var_ip_counts = args[2]
    all_counts = args[3]
    var_effi_x0 = args[4]
    #var_input_counts=all_counts-var_ip_counts

    sum_0=0.0; sum_1=0.0; sum_2=0.0

    for i in range(len(var_m6a_psi)):
        sum_0+=-1*(var_ip_counts[i]*var_group[i]*var_effi_x0-(all_counts[i]*var_group[i]*var_effi_x0*np.exp(beta_1*var_group[i]*var_effi_x0+\
                   beta_2*var_m6a_psi[i]+theta_0))/(1+np.exp(beta_1*var_group[i]*var_effi_x0+\
                   beta_2*var_m6a_psi[i]+theta_0)))

        sum_1+=-1*(var_ip_counts[i]*var_m6a_psi[i]-(all_counts[i]*var_m6a_psi[i]*np.exp(beta_1*var_group[i]*var_effi_x0+\
                   beta_2*var_m6a_psi[i]+theta_0))/(1+np.exp(beta_1*var_group[i]*var_effi_x0+\
                   beta_2*var_m6a_psi[i]+theta_0)))

        sum_2=-1*(var_ip_counts[i]-(all_counts[i]*np.exp(beta_1*var_group[i]*var_effi_x0+\
                   beta_2*var_m6a_psi[i]+theta_0))/(1+np.exp(beta_1*var_group[i]*var_effi_x0+\
                   beta_2*var_m6a_psi[i]+theta_0)))

    value_der=np.asarray([sum_0, sum_1, sum_2])

    return value_der

#############################################################################
##  likelihood with no parameter-contraint
#############################################################################
def MLE_no_constrain(beta_theta, group_id, psi, ip_counts, all_counts, effi_1):
    xopt=fmin_l_bfgs_b(likelihood_function, beta_theta, fprime=likelihood_fun_der,
                       args=[group_id, psi, ip_counts, all_counts, effi_1],
                       bounds=[[-0.05, 0.05], [None, None], [None, None]],
                       m=10, factr=1e5, iprint=-1, maxls=50)
    sum_value_0= likelihood_function(xopt[0], *[group_id, psi, ip_counts, all_counts, effi_1])
    return sum_value_0

#############################################################################
##  likelihood with parameter-contraint
#############################################################################
def MLE_constrain(beta_theta, group_id, psi, ip_counts, all_counts, effi_1):
    xopt=fmin_l_bfgs_b(likelihood_function, beta_theta, fprime=likelihood_fun_der,
                       args=[group_id, psi, ip_counts, all_counts, effi_1],
                       bounds=[[0, 0], [None, None], [None, None]],
                       m=10, factr=1e5, iprint=-1, maxls=50)
    sum_value_1= likelihood_function(xopt[0], *[group_id, psi, ip_counts, all_counts, effi_1])
    return sum_value_1

#############################################################################
##  pavalue adjustion for multiple test; BH-adjusted FDR were calculated
#############################################################################
def fdr_adjust(ls_pv):
    m = len(ls_pv)
    if m<1:
        return []
    pv=np.asfarray(ls_pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    qvalues = m * [0]
    mincoeff = pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues

#############################################################################
##  logit psi transformation
#############################################################################
def logit(IP, Input):
    ls_logit=[]
    if len(IP)==len(Input):
        for i in range(len(IP)):
            temp=(IP[i]+0.01)/(Input[i]+0.01)
            adds=(temp)/(1-temp+0.01)
            value=log(adds)
            ls_logit.append(value)
    return ls_logit


#############################################################################
##  function of likelihood ratio test; return the diff-test pvalue
#############################################################################
def likelihood_ratio_test(args):
    group_id = args[0]
    psi_ratio = args[1]
    ip_counts = args[2]
    input_counts=args[3]
    effi = args[4]

    yy = logit(ip_counts, input_counts)
    xx = []
    for i in range(len(psi_ratio)):
        xx.append([group_id[i], psi_ratio[i]])
    xx = np.asarray(xx)
    reg = linear_model.LinearRegression().fit(xx, yy)
    beta_theta = [reg.coef_[0], reg.coef_[1], 0] ## the paramaters from linear regression model were used as initial parameters

    MLE_value_0 = MLE_no_constrain(beta_theta, group_id, psi_ratio, ip_counts, input_counts, effi)
    MLE_value_1 = MLE_constrain(beta_theta, group_id, psi_ratio, ip_counts, input_counts, effi)

    pvalue = (1-scipy.stats.chi2.cdf(2*(MLE_value_1 - MLE_value_0), 1))

    if np.isnan(pvalue):
        pvalue = 1.0

    return pvalue


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='please provide the data file for differential m6A test!')
    parser.add_argument('--input', required=True, help='specify the input data file')
    parser.add_argument('--fdr', required=False, default=1.0, help='specify the fdr cutoff')
    parser.add_argument('--output', required=True, help='specify the output file')

    args = parser.parse_args()
    input_file=args.input
    fdr_cut=float(args.fdr)
    print (fdr_cut)
    output_file=args.output

    if len([input_file, output_file])<2:
        print("please specify your input and output file!")
        sys.exit(1)
    else:
        try:
           fread_data=open(input_file,'r')
           fread_data.readline()
           dic_temp={}
           ls_data_temp = []
        except EnvironmentError as err:
           raise NameError(str(err))

    for line in fread_data:
        temp = str(line).strip().split('\t')
        Peak_id = temp[0]
        ip_counts = np.asarray(str(temp[1]).strip().split(';'), dtype='int')
        input_counts = np.asarray(str(temp[2]).strip().split(';'), dtype='int')
        input_counts = ip_counts + input_counts

        psi=np.asarray(str(temp[3]).strip().split(';'), dtype='float')
        group_id=np.array(str(temp[4]).strip().split(';'), dtype='int')-1

        effi_1, effi_2, effi_3 = np.array(str(temp[5]).strip().split(';'), dtype='float')
        ls_data_temp.append([group_id, psi, ip_counts, input_counts, 1-effi_1])

        dic_temp_group={}
        for i in range(len(group_id)):
            dic_temp_group.setdefault(group_id[i], []).append(psi[i])

        group_key=sorted(dic_temp_group.keys())
        psi_1=np.mean(dic_temp_group[group_key[0]])
        psi_2=np.mean(dic_temp_group[group_key[1]])

        foldChange=log((psi_2+0.001)/(psi_1+0.001), 2)
        foldChange_1 = psi_2 - psi_1
        dic_temp.setdefault(Peak_id, [foldChange, foldChange_1])

    fread_data.close()

    pool = Pool()
    pvalue_ls = pool.map(likelihood_ratio_test, ls_data_temp)
    pool.close()
    pool.join()

    fdr = fdr_adjust(pvalue_ls)

    fread_data=open(input_file, 'r')
    fread_data.readline()
    fout_data=open(output_file, 'w')

    fout_data.writelines('Peak_id\tIP_counts\tInput_counts\tm6A_ratio\tSample_group_id\tBeta_level\tlogFC\tdiff_psi\tp-value\tFDR\n')
    i = 0
    for line in fread_data:
        temp = line.strip().split('\n')[0]
        key = str(temp).split('\t')[0]
        logFC = dic_temp[key][0]
        diff_pis = dic_temp[key][1]
        pvalue = pvalue_ls[i]

        if fdr[i]<= fdr_cut:
           # print fdr[i]
           fout_data.writelines('%s\t%.3f\t%.3f\t%.3e\t%.3e\n' % (temp, logFC, diff_pis, pvalue, fdr[i]))
        i += 1
    fread_data.close()
    fout_data.close()
