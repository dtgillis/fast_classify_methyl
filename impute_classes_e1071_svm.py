

import os
import numpy as np
from sklearn import cross_validation
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from scipy.stats import pearsonr
import argparse
import time
from rpy2.robjects import FactorVector
from rpy2.robjects import r


def clean_pearsonr(y_pred, y_test):

    pearson_r = pearsonr(y_pred, y_test)

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
    x_final = x
    y = np_y[:, 1]
    e1071 = importr('e1071')
    # make a tunning function using correlation
    error_func = r(
        'error.func <- function(x,y){ if (is.na(cor(as.numeric(x),as.numeric(y)))) return( 1.0 )  else return(1 - (cor(as.numeric(x),as.numeric(y)))^2) }')
    tune_func = e1071.tune_control(cross=5, error_fun=error_func)
    # start up R stuff
    print "starting r"
    print "fitting using cross validation"
    gamma = 2**np.arange(-1, 0, 1.0)
    cost = 2**np.arange(-1, 0, 1.0)
    # fit model using tune in R package.
    fit = e1071.tune_svm(x=x_final, y=FactorVector(y), gamma=gamma, cost=cost, tunecontrol=tune_func)
    #best_model = fit.rx2('best.model')
    total_time = time.time() - start_time
    #print "Finished fitting now cross validate for score"
    # now make a couple predictions using same data to get r squared
    # scores = []
    # for train_index, test_index in skf:
    #
    #     y_test = y[test_index]
    #     x_test = x_final[test_index]
    #
    #     y_pred = np.array(e1071.predict_svm(best_model, x_test))
    #
    #     scores.append(clean_pearsonr(y_pred, y_test))
    # what were the results ?
    best_parameters = fit.rx('best.parameters')[0]
    best_performance = 1.0 - float(fit.rx('best.performance')[0][0])
    gamma = float(best_parameters.rx('gamma')[0][0])
    cost = float(best_parameters.rx('cost')[0][0])
    # write out results.
    out = '{0:s} {1:f} {2:f} gamma {3:f} cost {4:f}'.format(args.snp_name[0], best_performance, total_time, gamma, cost)
    out += os.linesep
    f = open(args.out_file[0], mode='ab')
    f.write(out)
    f.close()


if __name__ == '__main__':
    main_program()






