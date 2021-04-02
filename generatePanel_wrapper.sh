#!/usr/bin/env bash

#! RUN : bash generatePanel_wrapper.sh <TAG> <CHR>

TAG=$1
CHR=$2

[[ -z "$TAG" ]] && { echo "ERROR: No TAG given for this run"; exit 1; }
[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }

source ${TAG}.config

uid=`date | md5sum | cut -c1-8`
RUN_NAME="${uid}"
mkdir -p "$RUN_NAME/logs"
cd $RUN_NAME
cp ../${TAG}.config .


# prepare files
# Take VCF value and look to see if need to make CHR specific VCF file

if [[ ! -e $VCF_ALL && -e $VCF ]]
then 
  echo "ERROR - Unable to find VCF file for ${TAG}";
  exit 1;
fi

# slurm_vcfChr.sh
jid1=$(sbatch -J ${TAG}.vcfChr ${HOME}/scripts/Imputation/slurm/slurm_vcfChr.sh ${TAG} ${CHR})

echo $jid1