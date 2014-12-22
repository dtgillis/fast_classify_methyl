import os
import numpy as np
import sys
from sklearn.grid_search import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import scale
from scipy.stats import pearsonr
import argparse
from sklearn.metrics import make_scorer
from sklearn.cross_validation import StratifiedKFold
import time
from scipy.io import savemat


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

    # make sure the directory exists
    if not os.path.exists(args.work_dir[0]):
        sys.exit(1)
    else:
        os.chdir(args.work_dir[0])
    start_time = time.time()
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
    #x = scale(x)
    x_final = x
    y = np_y[:, 1]
    a_first_split = .001
    a_range = 2**np.arange(-10.0, -2.0)
    print a_range
    loss = ['perceptron', 'hinge']
    l1_first_split = .05
    l1_range = np.arange(0, 1, l1_first_split)
    param_grid = dict(alpha=a_range, l1_ratio=l1_range, loss=loss)
    cv = StratifiedKFold(y=y, n_folds=5)
    clf = GridSearchCV(SGDClassifier(penalty='elasticnet', class_weight='auto'), param_grid=param_grid, cv=cv, scoring=make_scorer(clean_pearsonr), refit=True, n_jobs=7)
    clf.fit(X=x, y=y)
    end_time = time.time()

    total_time = end_time - start_time
    f = open(args.out_file[0], mode='ab')
    out = '{0:s} {1:f} {2:f} l1_ratio {3:f} alpha {4:f} loss {5:s}\n'.format(args.snp_name[0], clf.best_score_,
                                                                             total_time,
                                                                             clf.best_params_['l1_ratio'],
                                                                             clf.best_params_['alpha'],
                                                                             clf.best_params_['loss'])
    f.write(out)
    f.close()

if __name__ == '__main__':
    main_program()



