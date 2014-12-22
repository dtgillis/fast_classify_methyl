#!/bin/bash


RUN_DIR=
OUT_FILE=#out#
SNP_LIST=#snp#
#chrm number only
CHRM=#chrm#
PYTHON_SCRIPT=haladfpoj
BETA_FILE_DIR=/projects/sequence_analysis/vol2/dtgillis/methyl_impute
MAPPING_FILE=/projects/sequence_analysis/vol2/dtgillis/methyl_impute/Heart_Go/overlap.GWAS+methy.txt

while read snp

do

    # run plink to extract snp in 0,1,2 format

    plink --noweb --no-pheno --map3 --ped chr${CHRM}.extract.ped \
    --map chr${CHRM}.map  --recodeA --snp ${snp} --out ${snp}


    PLINK_FILE=${RUN_DIR}/${snp}.raw
    BETA_FILE=${BETA_FILE_DIR}/methyl.beta.chrm.$CHRM

    python $PYTHON_SCRIPT  --work_dir /home/dtgillis/winter_rotation/scratch_work \
    --snp_name $snp --out_file $OUT_FILE \
    --mapping_file $MAPPING_FILE \
    --beta_file $BETA_FILE \
    --plink_file $PLINK_FILE


done < $SNP_LIST