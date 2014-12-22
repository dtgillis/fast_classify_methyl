#!/bin/bash

n_jobs=10
MAIN_WORK=$PWD
SNP_LIST_DIR=/projects/sequence_analysis/vol2/dtgillis/methyl_impute/run_methyl/snp_lists

for chrm_num in 9
do
    if [ ! -d "$MAIN_WORK/runs/$chrm_num" ];
    then
    mkdir -p ${MAIN_WORK}/runs/$chrm_num
    mkdir -p ${MAIN_WORK}/runs/$chrm_num/submit_scripts
    fi
    #split the snp files into chunks to compute seperately
    split -d -n n_jobs ${SNP_LIST_DIR}/chrm${chrm_num}.snps.gemes ${SNP_LIST_DIR}/chrm${chrm_num}.split.
    RUN_SNP_LIST=${MAIN_WORK}/runs/${chrm_num}/snp_lists
    mkdir ${RUN_SNP_LIST}
    mv ${SNP_LIST_DIR}/chrm${chrm_num}.split.* ${RUN_SNP_LIST}/

    # make run file for each snp subset
    for snp_file in ${RUN_SNP_LIST}/chrm${chrm_num}.split.*
    do
        #make the run files
        sed -e "s%#snp_list#%${snp_file}%g" -e "s%#out_file#%chrm${chrm_num}.results.dat%g" \
        -e "s%#chrm#%${chrm_num}%g" run_classifier.sh \
        > ${MAIN_WORK}/runs/${chrm_num}/submit_scripts/${snp_file}.run.sh
        # make run file executable
        chmod +x ${MAIN_WORK}/runs/${chrm_num}/submit_scripts/${snp_file}.run.sh
        # setup a submit script for this run
        echo -e "qsub -q serial -o ${snp_file}.o -e ${snp_file}.e ${MAIN_WORK}/runs/${chrm_num}/submit_scripts/./${snp_file}.run.sh" \
        >> ${MAIN_WORK}/runs/${chrm_num}/submit_scripts/submit.sh
        echo -e "rm ${snp_file}.o ${snp_file}.e" >> ${MAIN_WORK}/runs/${chrm_num}/submit_scripts/submit.sh
    done;

    # make the submit script executable
    chmod +x ${MAIN_WORK}/runs/${chrm_num}/submit_scripts/submit.sh

done;
