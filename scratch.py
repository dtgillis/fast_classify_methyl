

import os
import numpy as np
import sys
import sklearn.linear_model as lin_model
from sklearn import cross_validation
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn import tree
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.preprocessing import normalize, scale
from sklearn.metrics.pairwise import cosine_similarity
from sklearn import neighbors
from sklearn.svm import l1_min_c
from sklearn.naive_bayes import GaussianNB
from sklearn.grid_search import GridSearchCV
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from sklearn.tree import DecisionTreeClassifier
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestClassifier


os.chdir('/home/dtgillis/winter_rotation/scratch_work')

beta_file = open('methyl.beta.chrm.9', 'r')
# get methyl sample names
samples = [sample.strip("\"") for sample in beta_file.readline().strip().split()]
total_samples = len(samples)

methyl_mat_lookup_dict = dict()
for i in range(0, total_samples):
    methyl_mat_lookup_dict[samples[i]] = i

methyl_to_gwas_tuples = [line.strip().split() for line in open('overlap.GWAS+methy.txt', 'r')]

gwas_methyl_lookup_dict = dict()

for lookup in methyl_to_gwas_tuples:
    gwas_methyl_lookup_dict[lookup[1]] = lookup[0]

beta_chrm_9 = np.genfromtxt('methyl.beta.chrm.9', skip_header=1, usecols=range(1, total_samples+1))

genotype_tuples = [(line.strip().split()[0], line.strip().split()[-1]) for line in open('plink.raw', 'r')]
# get rid of first header line
genotype_tuples = genotype_tuples[1:]
genotype_dict = dict()
x_class_0 = []
x_class_1 = []
x_class_2 = []
for genotype_tuple in genotype_tuples:

    if genotype_tuple[1] not in ["NA"]:
        # really awful dictionary stuff.
        if genotype_tuple[0] in gwas_methyl_lookup_dict:
            if gwas_methyl_lookup_dict[genotype_tuple[0]] in methyl_mat_lookup_dict:
                genotype_dict[methyl_mat_lookup_dict[
                    gwas_methyl_lookup_dict[genotype_tuple[0]]]] = int(genotype_tuple[1])


print len(genotype_dict)

valid_cols = [key for key in genotype_dict]
valid_cols = sorted(valid_cols)
sorted_y = sorted(genotype_dict.items(), key=lambda genotype: genotype[0])
np_y = np.array(sorted_y)

X = beta_chrm_9[:, valid_cols]
# similarity = cosine_similarity(X)
# plt.matshow(similarity)
# plt.colorbar()
# plt.show()

X = X.T
M = np.log2(np.divide(X, 1-X))
X = scale(X)
# similarity = cosine_similarity(X)
# #print similarity
# plt.matshow(similarity)
# plt.colorbar()
# plt.show()
# correlation = np.corrcoef(X)
# plt.matshow(correlation)
# plt.colorbar()
# plt.show()
# M = np.log2(np.divide(X, 1-X))
M = scale(M)
#X_final = np.ones((X.shape[0], X.shape[1]+1))
#M_final = np.ones((M.shape[0], M.shape[1] + 1))
M_final = M
X_final = X
y = np_y[:, 1]
skf = cross_validation.StratifiedKFold(y, n_folds=5)
score_beta_log = 0
total = 0
score_m = 0
score_svm = 0
score_knn = 0
score_knn_m = 0
score_m_log = 0
score_svc_m = 0
score_b_nb = 0
score_m_nb = 0
score_dtc_m = 0
score_dtc_b = 0
# smac = importr('smac')
# print smac
# C_range = 10.0 ** np.arange(-2, 9)
# gamma_range = 10.0 ** np.arange(-5, 4)
# param_grid = dict(gamma=gamma_range, C=C_range)
# cv = cross_validation.StratifiedKFold(y=y, n_folds=3)
# grid = GridSearchCV(SVC(), param_grid=param_grid, cv=cv)
# grid.fit(X, y)

# print("The best classifier is: ", grid.best_estimator_)

for train_index, test_index in skf:



    #get the training and test data for betas
    x_train, x_test = X_final[train_index], X_final[test_index],
    y_train, y_test = y[train_index], y[test_index]

    # get the training data for the M values
    x_train_m, x_test_m = M_final[train_index], M_final[test_index]

    # get min c values for l1 penalties
    # min_beta_c = l1_min_c(x_train, y_train)
    # min_m_c = l1_min_c(x_train_m, y_train)

    # print "starting r"
    # fit = smac.smac(x_train, y_train)
    # predict = smac.predict_smac(fit, x_test)
    # print predict

    # decision tree classifiers
    #beta values
    dtc = DecisionTreeClassifier()
    dtc.fit(x_train, y_train)
    # m values
    dtc_m = DecisionTreeClassifier()
    dtc_m.fit(x_train_m, y_train)

    # random forest classifier

    #beta values
    rfc = RandomForestClassifier(n_jobs=4)
    rfc.fit(x_train, y_train)

    # m values
    rfc_m = RandomForestClassifier(n_jobs=4)
    rfc_m.fit(x_train_m, y_train)

    #predict
    y_predict_rfc = rfc.predict(x_test)
    y_predict_rfc_m = rfc_m.predict(x_test_m)


    # predict dtcs
    y_predict_dtc_b = dtc.predict(x_test)
    y_predict_dtc_m = dtc_m.predict(x_test_m)

    # score dtcs
    score_dtc_b += metrics.accuracy_score(y_test, y_predict_dtc_b, normalize=False)
    score_dtc_m += metrics.accuracy_score(y_test, y_predict_dtc_m, normalize=False)

    # KNN
    # beta values
    knn = neighbors.KNeighborsClassifier(10, weights='distance')
    knn.fit(x_train, y_train)

    # m values
    knn_m = neighbors.KNeighborsClassifier(10, weights='distance')
    knn_m.fit(x_train_m, y_train)

    # predict knn classes
    y_predict_knn_m = knn_m.predict(x_test_m)
    y_predict_knn = knn.predict(x_test)
    score_knn += metrics.accuracy_score(y_test, y_predict_knn, normalize=False)
    score_knn_m += metrics.accuracy_score(y_test, y_predict_knn_m, normalize=False)

    # svc

    # beta values
    clf = SVC(C=1000)
    clf.fit(x_train, y_train)

    # m values
    clf_m = SVC(C=1000)
    clf_m.fit(x_train_m, y_train)

    # predict and score
    y_predict_svc_m = clf_m.predict(x_test_m)
    y_predict_svc = clf.predict(x_test)
    score_svc_m += metrics.accuracy_score(y_test, y_predict_svc_m, normalize=False)

    # Logistic Regression

    # beta values l2 penalty
    log_reg_beta_l2 = lin_model.LogisticRegression(C=1000, penalty='l1', class_weight='auto')
    log_reg_beta_l2.fit(x_train, y_train)

    # M
    log_reg_m_l2 = lin_model.LogisticRegression(C=1000, penalty='l1', class_weight='auto')
    log_reg_m_l2.fit(x_train_m, y_train)

    # score logistic models
    y_predict_log_m = log_reg_m_l2.predict(x_test_m)
    y_predict_beta = log_reg_beta_l2.predict(x_test)
    score_beta_log += metrics.accuracy_score(y_test, y_predict_beta, normalize=False)
    score_m_log += metrics.accuracy_score(y_test, y_predict_log_m, normalize=False)

    # rand log reg

    rand_reg_beta_l2 = lin_model.RandomizedLogisticRegression()
    rand_reg_beta_l2.fit(x_train, y_train)
    #print rand_reg_beta_l2.get_support()
    # naive bayes

    # beta values
    naive_beta = GaussianNB()
    naive_beta.fit(x_train, y_train)
    y_predict_nb_beta = naive_beta.predict(x_test)

    # m values
    naive_m = GaussianNB()
    naive_m.fit(x_train_m, y_train)
    y_predict_nb_m = naive_m.predict(x_test_m)

    score_m_nb += metrics.accuracy_score(y_test, y_predict_nb_m, normalize=False)
    score_b_nb += metrics.accuracy_score(y_test, y_predict_nb_beta, normalize=False)

    total += y_test.shape[0]

    score_svm += metrics.accuracy_score(y_test, y_predict_svc, normalize=False)
    print "Logistic report: ", metrics.classification_report(y_test, y_predict_beta)
    print "SVM report: ", metrics.classification_report(y_test, y_predict_svc)
    print "KNN report: ", metrics.classification_report(y_test, y_predict_knn)
    print "Logistic M: ", metrics.classification_report(y_test, y_predict_log_m)
    print "SVM M: ", metrics.classification_report(y_test, y_predict_svc_m)
    print "KNN M: ", metrics.classification_report(y_test, y_predict_knn_m)
    print "NB M: ", metrics.classification_report(y_test, y_predict_nb_m)
    print "NB B: ", metrics.classification_report(y_test, y_predict_nb_beta)
    print "DTC B: ", metrics.classification_report(y_test, y_predict_dtc_b)
    print "DTC M: ", metrics.classification_report(y_test, y_predict_dtc_m)
    print "Logistic r2: ", pearsonr(y_test, y_predict_beta)
    print "SVM r2: ", pearsonr(y_test, y_predict_svc)
    print "KNN r2: ", pearsonr(y_test, y_predict_knn)
    print "Logistic r2 M: ", pearsonr(y_test, y_predict_log_m)
    print "SVM M r2: ", pearsonr(y_test, y_predict_svc_m)
    print "KNN M r2: ", pearsonr(y_test, y_predict_knn_m)
    print "NB M r2: ", pearsonr(y_test, y_predict_nb_m)
    print "NB B r2: ", pearsonr(y_test, y_predict_nb_beta)
    print "DTC B r2: ", pearsonr(y_test, y_predict_dtc_b)
    print "DTC M r2: ", pearsonr(y_test, y_predict_dtc_m)
    print "RFC B r2: ", pearsonr(y_test, y_predict_rfc)
    print "RFC M r2: ", pearsonr(y_test, y_predict_rfc_m)

print "Logistic: ", score_beta_log / float(total)
print "SVM: ", score_svm / float(total)
print "KNN: ", score_knn / float(total)
print "Logistic M: ", score_m_log / float(total)
print "SVM M: ", score_svc_m / float(total)
print "KNN M: ", score_knn_m / float(total)
print "NB beta: ", score_b_nb / float(total)
print "NB M: ", score_m_nb / float(total)






