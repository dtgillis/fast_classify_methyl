

import os
import numpy as np
from sklearn import cross_validation
from sklearn.preprocessing import normalize, scale
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from scipy.stats import pearsonr
import argparse
import time


def clean_pearsonr(y_true, y_pred):

    pearson_r = pearsonr(y_true, y_pred)
    if np.isnan(pearson_r[0]):
        return .0
    else:
        return pearson_r[0]**2


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
    os.chdir(args.work_dir[0])
    start_time = time.time()
    beta_file = open(args.beta_file[0], 'r')
    # get methyl sample names
    samples = [sample.strip("\"") for sample in beta_file.readline().strip().split()]
    total_samples = len(samples)
    methyl_mat_lookup_dict = dict()
    for i in range(0, total_samples):
        methyl_mat_lookup_dict[samples[i]] = i

    methyl_to_gwas_tuples = [line.strip().split() for line in open(args.mapping_file[0], 'r')]

    gwas_methyl_lookup_dict = dict()

    for lookup in methyl_to_gwas_tuples:
        gwas_methyl_lookup_dict[lookup[1]] = lookup[0]

    beta_chrm_9 = np.genfromtxt(args.beta_file[0], skip_header=1, usecols=range(1, total_samples+1))

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
    sorted_y = sorted(genotype_dict.items(), key=lambda genotype: genotype[0])
    np_y = np.array(sorted_y)

    x = beta_chrm_9[:, valid_cols]
    x = x.T
    #m = np.log2(np.divide(x, 1-x))
    x = scale(x)
    #m = scale(m)
    #m_final = m
    x_final = x
    y = np_y[:, 1]
    skf = cross_validation.StratifiedKFold(y, n_folds=5)
    smac = importr('smac')
    print "starting r"
    print "fitting using cross validation"
    fit = smac.cvsmac(x_final, y)
    # predict = smac.predict_smac(fit, x_test)
    # predict_y = np.array(predict.rx('pred.y')).reshape((100, len(y_test)))
    # score_smac.append(clean_pearsonr(pearsonr(y_test, predict_y[-1])))
    min_error = fit.rx('min.error')[0][0]
    best_lambdas = fit.rx('best.lambda')[0]
    print fit.rx('best.lambda')
    total_time = time.time() - start_time
    pearson_r = []
    for train_index, test_index in skf:

        #get the training and test data for betas
        x_train, x_test = x_final[train_index], x_final[test_index]
        y_train, y_test = y[train_index], y[test_index]
        predictor = smac.predict_cvsmac(fit, x_test)
        best_y = predictor.rx('best.pred.y')[0][0]
        pearson_r.append(clean_pearsonr(y_test, best_y))

    out = '{0:s} {1:f} {2:f}'.format(args.snp_name[0], np.mean(pearson_r), total_time)
    for best_lambda in best_lambdas:
        out += ' {0:f}'.format(float(best_lambda))

    out += os.linesep
    f = open(args.out_file[0], mode='ab')
    f.write(out)
    f.close()


if __name__ == '__main__':
    main_program()






