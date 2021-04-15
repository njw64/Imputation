#!/usr/bin/env bash

#! RUN : bash generatePanel_wrapper.sh <CHR> <PROJECT>
#! Eg. : bash generatePanel_wrapper.sh 13 go_lab

CHR=$1
#PROJECT=$2

[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }
#[[ -z "$PROJECT" ]] && { echo "ERROR: No PROJECT provided for this run"; exit 1; }

JOBS=""
mkdir -p chr${CHR}/logs; cd chr${CHR};
cp ../setup.config .
source setup.config

mkdir $GWAS_PLINK; cd $GWAS_PLINK

jid1=$(sbatch -J ${GWAS_PLINK}.prep ${HOME}/scripts/Imputation/slurm/prepareGWAS.sh ${CHR})
jid2=$(sbatch -J ${GWAS_PLINK}.shapeit --dependency=afterok:${jid1##* } ${HOME}/scripts/Imputation/slurm/shapeitGWAS.sh ${CHR})
