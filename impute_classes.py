import os
import numpy as np
import sys
import sklearn.linear_model as lin_model
from sklearn import cross_validation
from sklearn.svm import SVC
from sklearn.preprocessing import scale
from sklearn import neighbors
from sklearn.naive_bayes import GaussianNB
from sklearn.tree import DecisionTreeClassifier
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestClassifier
import argparse
import io

def clean_pearsonr(pearson_r):

    if np.isnan(pearson_r[0]):
        return .0
    else:
        return pearson_r[0]


def main_program():

    parser = argparse.ArgumentParser(description='Impute some genotype classes')
    parser.add_argument('--work_dir', type=str, nargs=1, help='working directory')
    parser.add_argument('--snp_name', type=str, nargs=1, help='current snp being imputed')
    parser.add_argument('--out_file', type=str, nargs=1, help='name of file to write results to ')
    parser.add_argument('--beta_file', type=str, nargs=1, help='name of the file with beta values')
    parser.add_argument('--plink_file', type=str, nargs=1, help='Name of plink file coded 0,1,2')
    parser.add_argument('--mapping_file', type=str, nargs=1, help='File with mapping between sample names in methylation '
                                                                  'array data to genotype plink sample names. Methylation '
                                                                  'data sample names should come first')

    args = parser.parse_args()

    # make sure the directory exists
    if not os.path.exists(args.work_dir[0]):
        sys.exit(1)
    else:
        os.chdir(args.work_dir[0])

    beta_file = open(args.beta_file[0], 'r')
    # get methyl sample names
    samples = [sample.strip("\"") for sample in beta_file.readline().strip().split()]
    total_samples = len(samples)

    # this is mapping from methyl sample name to gwas sample name
    methyl_mat_lookup_dict = dict()
    for i in range(0, total_samples):
        methyl_mat_lookup_dict[samples[i]] = i

    methyl_to_gwas_tuples = [line.strip().split() for line in open(args.mapping_file[0], 'r')]

    gwas_methyl_lookup_dict = dict()

    for lookup in methyl_to_gwas_tuples:
        gwas_methyl_lookup_dict[lookup[1]] = lookup[0]
    #header is methyl sample names, and 0th column in the probe names
    beta_values = np.genfromtxt(args.beta_file[0], skip_header=1, usecols=range(1, total_samples+1))
    #genotype tuple is ( sample name (str) , genotype, (0,1,2))
    genotype_tuples = [(line.strip().split()[0], line.strip().split()[-1]) for line in open(args.plink_file[0], 'r')]
    # get rid of first header line
    genotype_tuples = genotype_tuples[1:]
    genotype_dict = dict()
    for genotype_tuple in genotype_tuples:

        if genotype_tuple[1] not in ["NA"]:
            # really awful dictionary stuff.
            if genotype_tuple[0] in gwas_methyl_lookup_dict:
                if gwas_methyl_lookup_dict[genotype_tuple[0]] in methyl_mat_lookup_dict:
                    genotype_dict[methyl_mat_lookup_dict[
                        gwas_methyl_lookup_dict[genotype_tuple[0]]]] = int(genotype_tuple[1])


    valid_cols = [key for key in genotype_dict]
    valid_cols = sorted(valid_cols)
    if len(valid_cols) == 0:
        print 'No valid samples with genotypes snp {0:s}'.format(args.snp_name[0])
        f = open(args.out_file[0], mode='a')
        out = '{1:s},{0:1.2f},{0:1.2f},{0:1.2f},{0:1.2f},{0:1.2f},{0:1.2f},' \
            '{0:1.2f},{0:1.2f},{0:1.2f},{0:1.2f},{0:1.2f},{0:1.2f}\n'.format(.0, args.snp_name[0])
        f.write(out)
        sys.exit(1)
    #sort by row number
    sorted_y = sorted(genotype_dict.items(), key=lambda genotype: genotype[0])
    np_y = np.array(sorted_y)

    # get the valid samples from beta file
    x = beta_values[:, valid_cols]
    x = x.T
    m = np.log2(np.divide(x, 1-x))
    x = scale(x)
    m = scale(m)
    m_final = m
    x_final = x
    y = np_y[:, 1]
    skf = cross_validation.StratifiedKFold(y, n_folds=5)
    score_log_b = []
    score_svm = []
    score_knn = []
    score_knn_m = []
    score_m_log = []
    score_svc_m = []
    score_b_nb = []
    score_m_nb = []
    score_dtc_m = []
    score_dtc_b = []
    score_rfc_b = []
    score_rfc_m = []
    for train_index, test_index in skf:
        #get the training and test data for betas
        x_train, x_test = x_final[train_index], x_final[test_index],
        y_train, y_test = y[train_index], y[test_index]

        # get the training data for the M values
        x_train_m, x_test_m = m_final[train_index], m_final[test_index]

        # decision tree classifiers
        #beta values
        dtc = DecisionTreeClassifier()
        dtc.fit(x_train, y_train)
        # m values
        dtc_m = DecisionTreeClassifier()
        dtc_m.fit(x_train_m, y_train)

        # predict dtcs
        y_predict_dtc_b = dtc.predict(x_test)
        y_predict_dtc_m = dtc_m.predict(x_test_m)


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

        # svc

        # beta values
        clf = SVC(C=1)
        clf.fit(x_train, y_train)

        # m values
        clf_m = SVC(C=1)
        clf_m.fit(x_train_m, y_train)

        # predict and score
        y_predict_svc_m = clf_m.predict(x_test_m)
        y_predict_svc = clf.predict(x_test)

        # Logistic Regression

        # beta values l2 penalty
        log_reg_beta_l2 = lin_model.LogisticRegression(C=1, penalty='l2', class_weight='auto')
        log_reg_beta_l2.fit(x_train, y_train)

        # M
        log_reg_m_l2 = lin_model.LogisticRegression(C=1, penalty='l2', class_weight='auto')
        log_reg_m_l2.fit(x_train_m, y_train)

        # score logistic models
        y_predict_log_m = log_reg_m_l2.predict(x_test_m)
        y_predict_beta = log_reg_beta_l2.predict(x_test)

        # naive bayes

        # beta values
        naive_beta = GaussianNB()
        naive_beta.fit(x_train, y_train)
        y_predict_nb_beta = naive_beta.predict(x_test)

        # m values
        naive_m = GaussianNB()
        naive_m.fit(x_train_m, y_train)
        y_predict_nb_m = naive_m.predict(x_test_m)

        score_log_b.append(clean_pearsonr(pearsonr(y_test, y_predict_beta)))
        score_svm.append(clean_pearsonr(pearsonr(y_test, y_predict_svc)))
        score_knn.append(clean_pearsonr(pearsonr(y_test, y_predict_knn)))
        score_m_log.append(clean_pearsonr(pearsonr(y_test, y_predict_log_m)))
        score_svc_m.append(clean_pearsonr(pearsonr(y_test, y_predict_svc_m)))
        score_knn_m.append(clean_pearsonr(pearsonr(y_test, y_predict_knn_m)))
        score_m_nb.append(clean_pearsonr(pearsonr(y_test, y_predict_nb_m)))
        score_b_nb.append(clean_pearsonr(pearsonr(y_test, y_predict_nb_beta)))
        score_dtc_b.append(clean_pearsonr(pearsonr(y_test, y_predict_dtc_b)))
        score_dtc_m.append(clean_pearsonr(pearsonr(y_test, y_predict_dtc_m)))
        score_rfc_b.append(clean_pearsonr(pearsonr(y_test, y_predict_rfc)))
        score_rfc_m.append(clean_pearsonr(pearsonr(y_test, y_predict_rfc_m)))

    f = open(args.out_file[0], mode='ab')
    out = '{12:s},{0:1.2f},{1:1.2f},{2:1.2f},{3:1.2f},{4:1.2f},{5:1.2f},' \
        '{6:1.2f},{7:1.2f},{8:1.2f},{9:1.2f},{10:1.2f},{11:1.2f}\n'.format(
        np.mean(score_log_b), np.mean(score_m_log), np.mean(score_svm), np.mean(score_svc_m),
        np.mean(score_knn), np.mean(score_knn_m), np.mean(score_rfc_b), np.mean(score_rfc_m),
        np.mean(score_dtc_b), np.mean(score_dtc_m), np.mean(score_b_nb), np.mean(score_m_nb),
        args.snp_name[0])
    f.write(out)
    f.close()

if __name__ == '__main__':
    main_program()



