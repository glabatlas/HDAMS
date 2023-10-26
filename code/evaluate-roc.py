import pandas as pd
import numpy as np
import os
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import mutual_info_classif
from sklearn.externals import joblib
from sklearn.metrics import roc_curve, auc, f1_score, precision_recall_curve
from sklearn.preprocessing import label_binarize
import warnings
# import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

def leave_one_out(fold_num, X_matrix, y_target, *classifiers):
    accuracy_score = []
    F1_score = []
    dic_micro = {}

    y_target_binary = label_binarize(y_target, classes=[1, 2, 3])
    n_classes = y_target_binary.shape[1]
    # print (y_target_binary)
    kf = StratifiedKFold(n_splits=fold_num, shuffle=False)
    mean_fpr = np.linspace(0, 1, 24)
    mean_precision =np.linspace(0, 1, 24)

    for train_index, test_index in kf.split(X_matrix, y_target):
        X_train, X_test = X_matrix[train_index], X_matrix[test_index]
        y_train, y_test = y_target[train_index], y_target[test_index]
        y_train_binary, y_test_binary = y_target_binary[train_index], y_target_binary[test_index]

        if len(classifiers)==1:
            classifier = classifiers[0]
            classifier.fit(X_train, y_train.ravel())
            # print (classifier.classes_)
            score_i = classifier.score(X_test, y_test)
            accuracy_score.append(score_i)
            y_pred = classifier.predict(X_test)
            y_pred_prob = classifier.predict_proba(X_test)

        elif len(classifiers)>1:
            score_i=0
            y_pred_prob = np.zeros_like(y_test_binary)

            for classifier in classifiers:
                classifier.fit(X_train, y_train.ravel())
                score_i = score_i + classifier.score(X_test, y_test)
                y_pred_prob = y_pred_prob + classifier.predict_proba(X_test)

            score_i = score_i/len(classifiers)
            y_pred_prob = y_pred_prob/len(classifiers)
            accuracy_score.append(score_i)

            y_pred=[]
            for i in range(len(y_pred_prob)):
                max_prob = max(y_pred_prob[i, :])
                class_index = list(y_pred_prob[i, :]).index(max_prob)+1
                y_pred.append(class_index)

        f1_score_i = f1_score(y_test, y_pred, average='micro')
        F1_score.append(f1_score_i)

        dic_fpr = dict()
        dic_tpr = dict()
        dic_auc = dict()

        dic_precision = dict()
        dic_recall = dict()
        dic_auprc = dict()

       #####  Micro-strategy based ROC ###
        for i in range (n_classes):
            dic_fpr[i], dic_tpr[i], thresholds_roc= roc_curve(y_test_binary[:, i], y_pred_prob[:, i])
            dic_auc[i] = auc(dic_fpr[i], dic_tpr[i])

            dic_precision[i], dic_recall[i], thresholds_prc = precision_recall_curve(y_test_binary[:, i], y_pred_prob[:, i])
            dic_auprc[i] = auc(dic_recall[i], dic_precision[i])

        dic_fpr["micro"], dic_tpr["micro"], tholds_roc = roc_curve(y_test_binary.ravel(), y_pred_prob.ravel())
        dic_tpr["micro"] = np.interp(mean_fpr, dic_fpr["micro"], dic_tpr["micro"])
        dic_auc["micro"] = auc(mean_fpr, dic_tpr["micro"])

        dic_precision["micro"], dic_recall["micro"], tholds_prc = precision_recall_curve(y_test_binary.ravel(), y_pred_prob.ravel())
        # aver_precision=average_precision_score(y_test_binary, y_pred_prob, average="micro")
        dic_recall["micro"] = np.interp(mean_precision, dic_precision["micro"], dic_recall["micro"])
        dic_auprc["micro"] = auc(dic_recall["micro"], mean_precision)

        # print(len(dic_recall["micro"]), len(dic_precision["micro"]))
        # print(dic_auprc)

        dic_micro.setdefault('fpr', []).append(dic_fpr['micro'])
        dic_micro.setdefault('tpr', []).append(dic_tpr['micro'])
        dic_micro.setdefault('auroc', []).append(dic_auc['micro'])

        dic_micro.setdefault('precision', []).append(dic_precision['micro'])
        dic_micro.setdefault('recall', []).append(dic_recall['micro'])
        dic_micro.setdefault('auprc', []).append(dic_auprc['micro'])


    acc_score_array = np.array(accuracy_score)
    F1_score_array = np.array(F1_score)

    micro_mean_tpr = np.mean(dic_micro['tpr'], axis=0)
    micro_mean_fpr = mean_fpr
    micro_mean_auc = auc(micro_mean_fpr, micro_mean_tpr)
    ls_micro_roc_res = (micro_mean_fpr, micro_mean_tpr, micro_mean_auc)
    # print(micro_mean_auc)
    # plt.plot(micro_mean_fpr, micro_mean_tpr)
    # plt.show()

    micro_mean_recall = np.mean(dic_micro['recall'], axis=0)
    micro_mean_precision = mean_precision
    micro_mean_auprc = auc(micro_mean_recall, micro_mean_precision)
    ls_micro_prc_res = (micro_mean_recall, micro_mean_precision, micro_mean_auprc)
    # print (micro_mean_auprc)
    # plt.plot(micro_mean_recall, micro_mean_precision)
    # plt.show()
    return acc_score_array, F1_score_array, ls_micro_roc_res, ls_micro_prc_res

def run_classifier_KFold(k, X_train, y_target):
    ftemp = open('test_roc_prc_new.txt', 'w')

    clf1 = SVC(C=1, kernel='linear', gamma='auto', probability=True)
    scores, scores_f1, ls_roc_res, ls_prc_res = leave_one_out(k, X_train, y_target, clf1)
    avg_accuracy = np.mean(scores)
    avg_f1 = np.mean(scores_f1)
    ls_fpr = ls_roc_res[0]
    ls_tpr = ls_roc_res[1]
    auroc = ls_roc_res[2]

    ls_recall= ls_prc_res[0]
    ls_precision= ls_prc_res[1]
    auprc = ls_prc_res[2]

    ftemp.writelines("frp\ttpr\tSVM\n")
    for fpr, tpr in zip(ls_fpr, ls_tpr):
        ftemp.writelines('%.3f\t%.3f\n'%(fpr, tpr))

    ftemp.writelines("\n")
    ftemp.writelines("recall\tprecision\tSVM\n")
    for reca, pres in zip(ls_recall, ls_precision):
       ftemp.writelines('%.3f\t%.3f\n'%(reca, pres))
    ftemp.writelines('SVM: acc->%.3f\tf1->%.3f\troc->%.3f\tprc->%.3f\n\n' % (avg_accuracy, avg_f1, auroc, auprc))

    print ('SVM: acc->%.3f; f1->%.3f; roc->%.3f; prc->%.3f' % (avg_accuracy, avg_f1, auroc, auprc))
    clf1_fit = clf1.fit(X_train, y_target)

    clf2 = LogisticRegression(C=10, penalty='l2', tol=1e-5, solver='lbfgs', multi_class='auto', dual=False, max_iter=500)
    scores, scores_f1, ls_roc_res, ls_prc_res = leave_one_out(k, X_train, y_target, clf2)
    avg_accuracy = np.mean(scores)
    avg_f1 = np.mean(scores_f1)
    ls_fpr = ls_roc_res[0]
    ls_tpr = ls_roc_res[1]
    auroc = ls_roc_res[2]

    ls_recall = ls_prc_res[0]
    ls_precision = ls_prc_res[1]
    auprc = ls_prc_res[2]

    ftemp.writelines("frp\ttpr\tLR\n")
    for fpr, tpr in zip(ls_fpr, ls_tpr):
        ftemp.writelines('%.3f\t%.3f\n' % (fpr, tpr))

    ftemp.writelines("\n")
    ftemp.writelines("recall\tprecision\tLR\n")
    for reca, pres in zip(ls_recall, ls_precision):
        ftemp.writelines('%.3f\t%.3f\n' % (reca, pres))
    ftemp.writelines("LR: acc->%.3f\tf1->%.3f\troc->%.3f\tprc->%.3f\n\n"%(avg_accuracy, avg_f1, auroc, auprc))

    print ('LR: acc->%.3f; f1->%.3f; roc->%.3f; prc->%.3f' % (avg_accuracy, avg_f1, auroc, auprc))
    clf2_fit=clf2.fit(X_train, y_target)

    scores, scores_f1, ls_roc_res, ls_prc_res = leave_one_out(k, X_train, y_target, clf1, clf2)
    avg_accuracy = np.mean(scores)
    avg_f1 =np.mean(scores_f1)

    ls_fpr = ls_roc_res[0]
    ls_tpr = ls_roc_res[1]
    auroc = ls_roc_res[2]

    ls_recall = ls_prc_res[0]
    ls_precision = ls_prc_res[1]
    auprc = ls_prc_res[2]

    ftemp.writelines("frp\ttpr\tSVM+LR\n")
    for fpr, tpr in zip(ls_fpr, ls_tpr):
        ftemp.writelines('%.3f\t%.3f\n' % (fpr, tpr))

    ftemp.writelines("\n")
    ftemp.writelines("recall\tprecision\tSVM+LR\n")
    for reca, pres in zip(ls_recall, ls_precision):
        ftemp.writelines('%.3f\t%.3f\n' % (reca, pres))
    ftemp.writelines("SVM+LR: acc->%.3f\tf1->%.3f\troc->%.3f\tprc->%.3f\n\n"%(avg_accuracy, avg_f1, auroc, auprc))
    ftemp.close()
    print('SVM+LR: acc->%.3f; f1->%.3f; roc->%.3f; prc->%.3f' % (avg_accuracy, avg_f1, auroc, auprc))

    return [clf1_fit, clf2_fit]

def classifier(data_feature, label_list, kfold):
    print ("supervised classifier based on features : \n")
    y_target = label_list
    X_train = np.array(data_feature)
    res_mode = run_classifier_KFold(k=kfold, X_train=X_train, y_target=y_target)
    return res_mode

def feature_selection(data_feature, data_lable):
    indices_seq_id = set()
    clf0 = ExtraTreesClassifier(n_estimators=101).fit(data_feature, data_lable)
    clf1 = LinearSVC(C=1, penalty="l1", dual=False, max_iter=20000).fit(data_feature, data_lable)
    clf2 = mutual_info_classif(data_feature, data_lable, n_neighbors=11)

    model0 = SelectFromModel(clf0, prefit=True, threshold=0.002)
    indices0 = model0.get_support(indices=True)

    indices_seq_id = indices_seq_id.union((indices0))

    model1 = SelectFromModel(clf1, prefit=True)
    indices1 = model1.get_support(indices=True)

    indices_seq_id = indices_seq_id.union((indices1))

    indices2 = [i for i in range(len(clf2)) if clf2[i]>0.1]
    indices_seq_id = list(indices_seq_id.union((indices2)))

    return indices_seq_id


if __name__ =="__main__":
    if not os.path.exists("./model/"):
        os.makedirs("./model/")
    data = pd.read_csv('./temp/seq_feature_combined.txt', sep='\t', header = None)
    sample_id = data.iloc[:, 1]
    data_lable = np.asarray(data.iloc[:, 0])
    data_feature = np.asarray(data.iloc[:, 2:])
    fout_index = open('./training_data/selected_feature_index_test.txt', 'w')
    indices_seq_id = sorted(feature_selection(data_feature, data_lable))

    sample_info = pd.DataFrame(data.iloc[:, 0:2])
    selected_data = pd.DataFrame(data_feature[:, indices_seq_id])
    selected_data_sample = pd.concat([sample_info, selected_data.reset_index(drop=True)], axis=1)

    for x in indices_seq_id:
        fout_index.writelines('%s\t'% x)
    fout_index.writelines('\n')
    fout_index.close()
    np.savetxt('./training_data/training_data_selected_feature_test.txt', selected_data_sample, fmt='%s', delimiter='\t')

    fitting_model = classifier(selected_data, data_lable, 5)
    path = "./model/"
    joblib.dump(fitting_model[0], path + 'svm_test.joblib', compress=9)
    joblib.dump(fitting_model[1], path + 'lr_test.joblib', compress=9)
