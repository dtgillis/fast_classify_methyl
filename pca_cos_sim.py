

import os
import numpy as np
from sklearn import cross_validation
from sklearn.preprocessing import normalize, scale
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import PCA
from scipy.stats import pearsonr


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

genotype_tuples = [(line.strip().split()[0], line.strip().split()[-1]) for line in open('rs3203916.raw', 'r')]
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
for train_index, test_index in skf:



    #get the training and test data for betas
    x_train, x_test = X_final[train_index], X_final[test_index],
    y_train, y_test = y[train_index], y[test_index]

    # get the training data for the M values
    x_train_m, x_test_m = M_final[train_index], M_final[test_index]
    class_dict = dict()
    for genotype in [0, 1, 2]:
        extract = np.where(y_train == genotype)
        if len(extract[0]) > 0:
            class_dict[genotype] = x_train[extract[0]]

    mean_dict = dict()

    for genotype in class_dict:

        mean_x = np.mean(class_dict[genotype], axis=0)
        mean_dict[genotype] = mean_x

    cosine_sim = np.zeros((x_test.shape[0], 3))
    for i in range(x_test.shape[0]):

        tmp_x = x_test[i, :]
        for genotype, mean_x in mean_dict.items():
            cosine_sim[i, genotype] = cosine_similarity(tmp_x, mean_dict[genotype])



    y_predict = []
    for row in cosine_sim:

        y_predict.append(np.where(row == np.max(row))[0][0])

    print len(np.where((np.array(y_predict) - y_test) == 0)[0])/float(len(y_predict))
    print pearsonr(np.array(y_predict), y_test)














