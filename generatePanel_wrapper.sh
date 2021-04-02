#!/usr/bin/env bash

#! RUN : bash generatePanel_wrapper.sh <TAG> <CHR>

TAG=$1
CHR=$2

[[ -z "$TAG" ]] && { echo "ERROR: No TAG given for this run"; exit 1; }
[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }

source ${TAG}.config

mkdir -p "$ID/logs"; cd $ID
cp ../${TAG}.config .


# prepare files
# Take VCF value and look to see if need to make CHR specific VCF file

if [[ ! -e $VCF_ALL && -e $VCF ]]
then 
  echo "ERROR - Unable to find VCF file for ${TAG}";
  exit 1;
fi

# vcfChr.sh
jid1=$(sbatch -J ${TAG}.vcfChr ${HOME}/scripts/Imputation/slurm/vcfChr.sh ${TAG} ${CHR})
# vcfQC.sh - convert VCF to plink; perform QC; check SNPs; convert back to VCF
jid2=$(sbatch -J ${TAG}.vcfQC --dependency=afterok:${jid1##* } ${HOME}/scripts/Imputation/slurm/vcfQC.sh ${TAG} ${CHR})

# The construct ${jid1##* } isolates the last word
echo ${jid1##* }