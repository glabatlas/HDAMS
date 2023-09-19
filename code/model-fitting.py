import pandas as pd
import numpy as np
import os
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.svm import SVC, LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel, mutual_info_classif
from sklearn.model_selection import StratifiedKFold
from sklearn.externals import joblib

def leave_one_out(fold_num, X_matrix, y_target, classifier):
    #from sklearn.model_selection import KFold
    import sklearn.metrics as metrics
    accuracy_score = []
    F1_score = []
    kf = StratifiedKFold(n_splits=fold_num, shuffle=False)
    for train_index, test_index in kf.split(X_matrix, y_target):
        X_train, X_test = X_matrix[train_index], X_matrix[test_index]
        y_train, y_test = y_target[train_index], y_target[test_index]
        classifier.fit(X_train, y_train.ravel())
        score_i = classifier.score(X_test, y_test)
        accuracy_score.append(score_i)
        y_pred = classifier.predict(X_test)
        F1_score_i = metrics.f1_score(y_test, y_pred, average='macro')
        F1_score.append(F1_score_i)

    acc_score_matrix = np.array(accuracy_score)
    F1_score_matrix = np.array(F1_score)
    return acc_score_matrix, F1_score_matrix

def run_classifier_KFold(k, X_train, y_target):
    clf1 = SVC(C=1, kernel='linear', gamma='auto', probability=True)
    scores, scores_f1 = leave_one_out(k, X_train, y_target, clf1)
    avg_accuracy = np.mean(scores)
    avg_f1 = np.mean(scores_f1)
    print ('SVM: acc->%f; f1->%f' % (avg_accuracy, avg_f1))
    clf1_fit = clf1.fit(X_train, y_target)

    clf2 = LogisticRegression(C=10, penalty='l2', tol=1e-5, solver='lbfgs', multi_class='auto', dual=False, max_iter=500)
    scores, scores_f1 = leave_one_out(k, X_train, y_target, clf2)
    avg_accuracy = np.mean(scores)
    avg_f1 = np.mean(scores_f1)
    print ('LR: acc->%f; f1->%f' % (avg_accuracy, avg_f1))
    clf2_fit=clf2.fit(X_train, y_target)

    return [clf1_fit, clf2_fit]

def classifier(data_feature, label_list, kfold):
    print ("supervised classifier based on features : \n")
    y_target = label_list
    X_train = np.array(data_feature)
    res_mode = run_classifier_KFold(k=kfold, X_train=X_train, y_target=y_target)
    return res_mode

def feature_selection(data_feature, data_lable):
    indices_seq_id= set()
    clf0= ExtraTreesClassifier(n_estimators=101).fit(data_feature, data_lable)
    clf1= LinearSVC(C=1, penalty="l1", dual=False, max_iter=20000).fit(data_feature, data_lable)
    clf2= mutual_info_classif(data_feature, data_lable, n_neighbors=11)

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
    data = pd.read_csv('./temp/seq_feature_combined.txt', sep='\t', header=None)
    sample_id = data.iloc[:, 1]
    data_lable = np.asarray(data.iloc[:, 0])
    # print (len(data_lable))

    data_feature = np.asarray(data.iloc[:, 2:])
    fout_index = open('./training_data/selected_feature_index.txt', 'w')
    indices_seq_id = sorted(feature_selection(data_feature, data_lable))
    sample_info = pd.DataFrame(data.iloc[:, 0:2])

    selected_data = pd.DataFrame(data_feature[:, indices_seq_id])
    selected_data_sample = pd.concat([sample_info, selected_data.reset_index(drop=True)], axis=1)

    for x in indices_seq_id:
        fout_index.writelines('%s\t'% x)
    fout_index.writelines('\n')
    fout_index.close()

    np.savetxt('./training_data/training_data_selected_feature.txt', selected_data_sample, fmt='%s', delimiter='\t')
    fitting_model = classifier(selected_data, data_lable, 5)
    path = "./model/"
    joblib.dump(fitting_model[0], path + 'svm.joblib', compress=9)
    joblib.dump(fitting_model[1], path + 'lr.joblib', compress=9)
